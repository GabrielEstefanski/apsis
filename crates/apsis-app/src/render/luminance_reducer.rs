//! GPU-side luminance metering for auto-exposure.
//!
//! Reduces the HDR scene target to a single scalar — `mean(L^p)` across
//! every pixel, stored in the final 1×1 mip — which the CPU side then
//! takes the 1/p root of to obtain a soft-max luminance (see
//! [`crate::render::exposure`] for the rationale behind the power-mean
//! approach).
//!
//! # Pipeline
//!
//! Two render pipelines, both fullscreen-triangle fragment passes:
//!
//! 1. **init** — samples the HDR colour target, computes
//!    `L = dot(rgb, LUMA_WEIGHTS)`, writes `L^p` to mip 0 of a private
//!    `R16Float` reducer texture. Rendered at **half the HDR resolution**:
//!    the downstream meter doesn't need per-pixel precision and the 4×
//!    memory saving buys a level of the mip chain for free.
//!
//! 2. **reduce** — one draw per mip transition (mip N → mip N+1). Samples
//!    the previous mip with bilinear filtering; a single tap at the texel
//!    centre of the destination fetches the 4-pixel average of the source
//!    thanks to hardware filtering — same cost as four manual taps but
//!    one draw-call per level instead of four.
//!
//! The final 1×1 mip is copied into a small staging buffer and mapped
//! asynchronously; the CPU polls for completion between frames. Latency
//! is 1–2 frames, acceptable for an exposure EMA that smooths over ~0.5 s.
//!
//! # Resize tracking
//!
//! The mip chain depth depends on the HDR size. [`LuminanceReducer::ensure_size`]
//! rebuilds the texture and per-mip views when the scene size changes.
//! Pipelines and the staging buffer are size-independent.
//!
//! # Readback robustness
//!
//! Three staging buffers cycle so that a late-mapping buffer doesn't
//! block the next frame's submission. Each buffer tracks its own
//! mapping state via `Arc<Mutex<…>>` so the async callback can update
//! it without racing the main thread.

use std::mem::size_of;
use std::sync::{Arc, Mutex};

use bytemuck::{Pod, Zeroable};

use crate::render::exposure::{LUMA_WEIGHTS, SOFT_MAX_P};
use crate::render::hdr::HdrTarget;

/// Single-channel, half-precision float: enough dynamic range for
/// `L^p` values (L typically 0–10, p = 4 → L^p up to ~10⁴) without
/// the bandwidth of 32-bit.
const REDUCER_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::R16Float;

/// Number of staging buffers cycled for async readback. Three is the
/// classic pattern: frame N writes buffer i, frame N+1 writes (i+1),
/// frame N+2 reads back buffer i (now definitely unmapped). Two would
/// work in principle but leaves no slack if the GPU runs a frame
/// behind; four is wasteful.
const READBACK_RING: usize = 3;

/// One f16 value, padded to the minimum copy-buffer alignment (256 B
/// per row for `copy_texture_to_buffer`). The row itself is mostly
/// padding; we only read the first two bytes.
const READBACK_ROW_BYTES: u64 = 256;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct InitUniform {
    /// Rec.709 luma weights plus the p-norm exponent — packed together
    /// because they're the only per-frame-constant values the init
    /// shader needs. `_pad` brings the struct to 16 B for WGSL uniform
    /// alignment rules.
    luma_weights: [f32; 3],
    p_norm: f32,
}

/// One slot in the readback ring.
struct ReadbackSlot {
    buffer: wgpu::Buffer,
    /// Outer `Option<f32>` — `Some` once the async map callback has
    /// written a value; `None` means "mapping still in flight" or
    /// "never submitted". Wrapped in `Arc<Mutex<…>>` because the
    /// callback runs on whichever thread wgpu's async machinery
    /// services and must publish the value safely.
    pending: Arc<Mutex<Option<f32>>>,
    /// True between the `copy_texture_to_buffer` submission and the
    /// time the CPU reads the mapped bytes. Prevents us from asking
    /// wgpu to re-map an already-pending buffer.
    in_flight: bool,
}

impl ReadbackSlot {
    fn new(device: &wgpu::Device, label: &str) -> Self {
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size: READBACK_ROW_BYTES,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        Self { buffer, pending: Arc::new(Mutex::new(None)), in_flight: false }
    }
}

pub struct LuminanceReducer {
    // ── GPU pipelines ─────────────────────────────────────────────────────────
    init_pipeline: wgpu::RenderPipeline,
    reduce_pipeline: wgpu::RenderPipeline,
    sampler: wgpu::Sampler,

    /// Uniform buffer holding [`InitUniform`]. Re-uploaded each frame
    /// (constant values but cheap) — could be upload-once but keeping
    /// the tuning knobs live for future CPU-side adjustment is worth
    /// the 16 B write.
    init_uniform_buf: wgpu::Buffer,

    init_bgl: wgpu::BindGroupLayout,
    reduce_bgl: wgpu::BindGroupLayout,
    /// Bind group for the init pass; rebuilt when the HDR view changes.
    init_bg: Option<wgpu::BindGroup>,
    /// Last observed HDR generation. Triggers `init_bg` rebuild.
    init_bound_generation: u64,

    // ── Mip-chain texture ────────────────────────────────────────────────────
    mip_tex: Option<wgpu::Texture>,
    /// One view per mip level. `mip_views[0]` is the render target for
    /// `init`; `mip_views[i+1]` is the target for the reduce pass that
    /// samples `mip_views[i]`. We need separate views because a view
    /// used as an attachment cannot also be sampled in the same pass.
    mip_views: Vec<wgpu::TextureView>,
    /// Per-mip bind groups for the reduce pass. `reduce_bgs[i]` samples
    /// mip i and is consumed by the pass that renders into mip i+1.
    reduce_bgs: Vec<wgpu::BindGroup>,
    /// HDR-texture size the mip chain was sized for. When it changes
    /// we rebuild.
    current_hdr_size: [u32; 2],

    // ── Readback ring ────────────────────────────────────────────────────────
    readback: [ReadbackSlot; READBACK_RING],
    /// Index of the slot scheduled to receive the *next* copy. Advances
    /// modulo [`READBACK_RING`] after every submission.
    next_slot: usize,
    /// Last decoded `mean(L^p)` — carried forward when no new reading
    /// is available this frame. `None` until the first successful
    /// readback, which suppresses the first-frame exposure spike.
    last_reading: Option<f32>,
}

impl LuminanceReducer {
    pub fn new(device: &wgpu::Device) -> Self {
        let init_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("lum_reducer::init_bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        // Reduce pass has no uniform — binding set is just {texture, sampler}.
        let reduce_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("lum_reducer::reduce_bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
            ],
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("lum_reducer::shader"),
            source: wgpu::ShaderSource::Wgsl(REDUCER_SHADER.into()),
        });

        let init_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("lum_reducer::init_layout"),
            bind_group_layouts: &[Some(&init_bgl)],
            immediate_size: 0,
        });
        let reduce_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("lum_reducer::reduce_layout"),
            bind_group_layouts: &[Some(&reduce_bgl)],
            immediate_size: 0,
        });

        let init_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("lum_reducer::init_pipeline"),
            layout: Some(&init_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_fullscreen"),
                buffers: &[],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_init"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: REDUCER_FORMAT,
                    blend: None,
                    write_mask: wgpu::ColorWrites::RED,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                ..Default::default()
            },
            depth_stencil: None,
            multisample: Default::default(),
            multiview_mask: None,
            cache: None,
        });

        let reduce_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("lum_reducer::reduce_pipeline"),
            layout: Some(&reduce_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_fullscreen"),
                buffers: &[],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_reduce"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: REDUCER_FORMAT,
                    blend: None,
                    write_mask: wgpu::ColorWrites::RED,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                ..Default::default()
            },
            depth_stencil: None,
            multisample: Default::default(),
            multiview_mask: None,
            cache: None,
        });

        let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("lum_reducer::sampler"),
            address_mode_u: wgpu::AddressMode::ClampToEdge,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            address_mode_w: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            mipmap_filter: wgpu::MipmapFilterMode::Nearest,
            ..Default::default()
        });

        let init_uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("lum_reducer::init_uniform"),
            size: size_of::<InitUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let readback = [
            ReadbackSlot::new(device, "lum_reducer::readback_0"),
            ReadbackSlot::new(device, "lum_reducer::readback_1"),
            ReadbackSlot::new(device, "lum_reducer::readback_2"),
        ];

        Self {
            init_pipeline,
            reduce_pipeline,
            sampler,
            init_uniform_buf,
            init_bgl,
            reduce_bgl,
            init_bg: None,
            init_bound_generation: 0,
            mip_tex: None,
            mip_views: Vec::new(),
            reduce_bgs: Vec::new(),
            current_hdr_size: [0, 0],
            readback,
            next_slot: 0,
            last_reading: None,
        }
    }

    /// Ensures the mip-chain texture matches the current HDR size.
    /// Called at the start of each encode; cheap when the size hasn't
    /// changed.
    fn ensure_size(&mut self, device: &wgpu::Device, hdr_size: [u32; 2]) {
        if self.mip_tex.is_some() && self.current_hdr_size == hdr_size {
            return;
        }

        // Init renders at half resolution — plenty for metering and it
        // saves a mip level of reduce work. Clamp to ≥ 1 so small canvas
        // sizes don't break the chain.
        let w0 = (hdr_size[0] / 2).max(1);
        let h0 = (hdr_size[1] / 2).max(1);
        let max_dim = w0.max(h0);
        let mip_level_count = 32u32.saturating_sub(max_dim.leading_zeros()).max(1);

        let tex = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("lum_reducer::mip_chain"),
            size: wgpu::Extent3d { width: w0, height: h0, depth_or_array_layers: 1 },
            mip_level_count,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: REDUCER_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                | wgpu::TextureUsages::TEXTURE_BINDING
                | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });

        let mut views = Vec::with_capacity(mip_level_count as usize);
        for mip in 0..mip_level_count {
            views.push(tex.create_view(&wgpu::TextureViewDescriptor {
                label: Some("lum_reducer::mip_view"),
                base_mip_level: mip,
                mip_level_count: Some(1),
                ..Default::default()
            }));
        }

        // Build per-transition reduce bind groups: bg[i] samples views[i],
        // to be consumed by the pass that renders into views[i+1].
        let mut reduce_bgs = Vec::with_capacity(mip_level_count as usize - 1);
        for i in 0..(mip_level_count as usize - 1) {
            reduce_bgs.push(device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("lum_reducer::reduce_bg"),
                layout: &self.reduce_bgl,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: wgpu::BindingResource::TextureView(&views[i]),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: wgpu::BindingResource::Sampler(&self.sampler),
                    },
                ],
            }));
        }

        self.mip_tex = Some(tex);
        self.mip_views = views;
        self.reduce_bgs = reduce_bgs;
        self.current_hdr_size = hdr_size;
    }

    fn ensure_init_bg(&mut self, device: &wgpu::Device, hdr: &HdrTarget) {
        if self.init_bg.is_some() && self.init_bound_generation == hdr.generation() {
            return;
        }
        self.init_bg = Some(device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("lum_reducer::init_bg"),
            layout: &self.init_bgl,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    // Reflective plane only — luminous bypasses metering.
                    resource: wgpu::BindingResource::TextureView(hdr.view_r()),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::Sampler(&self.sampler),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.init_uniform_buf.as_entire_binding(),
                },
            ],
        }));
        self.init_bound_generation = hdr.generation();
    }

    /// Encodes the init + reduce-chain + 1×1-copy passes on `encoder`.
    /// Must be called **after** the scene has been drawn into `hdr` so
    /// the HDR texture's contents are ready to be sampled.
    pub fn encode(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        encoder: &mut wgpu::CommandEncoder,
        hdr: &HdrTarget,
    ) {
        self.ensure_size(device, hdr.size());
        self.ensure_init_bg(device, hdr);

        queue.write_buffer(
            &self.init_uniform_buf,
            0,
            bytemuck::bytes_of(&InitUniform { luma_weights: LUMA_WEIGHTS, p_norm: SOFT_MAX_P }),
        );

        let Some(init_bg) = self.init_bg.as_ref() else { return };

        // ── init pass: HDR → mip 0 of reducer (L^p) ──────────────────────────
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("lum_reducer::init_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.mip_views[0],
                    resolve_target: None,
                    depth_slice: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_pipeline(&self.init_pipeline);
            pass.set_bind_group(0, init_bg, &[]);
            pass.draw(0..3, 0..1);
        }

        // ── reduce chain ────────────────────────────────────────────────────
        for i in 0..self.reduce_bgs.len() {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("lum_reducer::reduce_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.mip_views[i + 1],
                    resolve_target: None,
                    depth_slice: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_pipeline(&self.reduce_pipeline);
            pass.set_bind_group(0, &self.reduce_bgs[i], &[]);
            pass.draw(0..3, 0..1);
        }

        // ── copy final 1×1 mip into the next free readback slot ─────────────
        let Some(tex) = self.mip_tex.as_ref() else { return };
        let last_mip = (self.mip_views.len() as u32).saturating_sub(1);
        let slot_idx = self.next_slot;
        let slot = &mut self.readback[slot_idx];
        // If this slot is still in-flight from two frames ago, the
        // previous map_async hasn't fired — skip the copy and try the
        // next slot on the next frame. In practice the ring is long
        // enough that this never trips outside of a paused/stuttering
        // application, but it's a safe fallback.
        if !slot.in_flight {
            encoder.copy_texture_to_buffer(
                wgpu::TexelCopyTextureInfo {
                    texture: tex,
                    mip_level: last_mip,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                wgpu::TexelCopyBufferInfo {
                    buffer: &slot.buffer,
                    layout: wgpu::TexelCopyBufferLayout {
                        offset: 0,
                        // `bytes_per_row = None` is allowed for single-row
                        // copies but some backends require an explicit value
                        // ≥ 256. We use the padded row to stay portable.
                        bytes_per_row: Some(READBACK_ROW_BYTES as u32),
                        rows_per_image: Some(1),
                    },
                },
                wgpu::Extent3d { width: 1, height: 1, depth_or_array_layers: 1 },
            );
            slot.in_flight = true;
            self.next_slot = (self.next_slot + 1) % READBACK_RING;
        }
    }

    /// Submits map requests for any in-flight slots and harvests any
    /// completed ones. Call once per frame **after** the frame's
    /// command buffer has been submitted to the queue (so the copy
    /// destination is actually populated).
    ///
    /// Returns the most recent successfully decoded `mean(L^p)` if
    /// new data arrived this call; otherwise `None`. The backend
    /// feeds this into [`crate::render::exposure::ExposureState::tick`]
    /// after converting it to a soft-max luminance via
    /// [`crate::render::exposure::decode_reduced_texel`].
    ///
    /// ## Per-slot state machine (stored in `pending`)
    ///
    /// * `None`                   — idle (never submitted, or just recycled)
    /// * `Some(f32::INFINITY)`    — `map_async` called, awaiting callback
    /// * `Some(f32::NAN)`         — callback fired OK, bytes ready to read
    ///
    /// The `INFINITY` sentinel must be written **before** calling
    /// `map_async`: without it, a second `poll()` call before the
    /// callback fires would see `pending == None` and try to map the
    /// same buffer again — which wgpu rejects with "already mapped".
    pub fn poll(&mut self, device: &wgpu::Device) -> Option<f32> {
        // 1) Request maps for any newly-in-flight slots that don't yet
        //    have a map request registered. The `INFINITY` sentinel
        //    latches the "map requested" state so we never re-call
        //    map_async on the same slot before its callback runs.
        for slot in self.readback.iter_mut() {
            if !slot.in_flight {
                continue;
            }
            let mut guard = slot.pending.lock().unwrap();
            if guard.is_none() {
                *guard = Some(f32::INFINITY);
                drop(guard);
                let pending = Arc::clone(&slot.pending);
                let buf_slice = slot.buffer.slice(..);
                buf_slice.map_async(wgpu::MapMode::Read, move |res| {
                    let mut guard = pending.lock().unwrap();
                    if res.is_ok() {
                        *guard = Some(f32::NAN);
                    } else {
                        // Map failed: clear so the slot can be recycled
                        // once the ring wraps. Publishing a bogus value
                        // would poison the EMA for several frames.
                        *guard = None;
                    }
                });
            }
        }

        // Drive wgpu's internal async machinery so callbacks have a
        // chance to fire. Poll (non-blocking) rather than Wait so we
        // never stall the frame — if the GPU hasn't finished yet, we'll
        // just pick the value up next call.
        let _ = device.poll(wgpu::PollType::Poll);

        // 2) For any slot whose callback has fired (sentinel = NaN),
        //    read the mapped bytes and record the real value.
        let mut newest: Option<f32> = None;
        for slot in self.readback.iter_mut() {
            if !slot.in_flight {
                continue;
            }
            let is_ready = matches!(*slot.pending.lock().unwrap(), Some(v) if v.is_nan());
            if !is_ready {
                continue;
            }

            // Read the first 2 bytes (f16) from the mapped buffer.
            let value = {
                let view = slot.buffer.slice(..).get_mapped_range();
                let bytes: [u8; 2] = [view[0], view[1]];
                f16_to_f32(u16::from_le_bytes(bytes))
            };
            slot.buffer.unmap();
            slot.in_flight = false;
            *slot.pending.lock().unwrap() = None;

            // Take the most recent reading if multiple slots matured in
            // the same frame (rare but possible on a stall).
            newest = Some(value);
        }

        if let Some(v) = newest {
            self.last_reading = Some(v);
        }
        newest
    }

    /// Most recent successfully decoded `mean(L^p)` — useful as a
    /// fallback when [`poll`](Self::poll) returns `None` this frame.
    #[inline]
    pub fn last_reading(&self) -> Option<f32> {
        self.last_reading
    }
}

/// Manual half-precision decode. `bytemuck` doesn't ship `f16` and we
/// only need one value per frame, so a 20-line scalar decoder is
/// lighter than pulling a crate. Based on the IEEE 754 binary16 layout:
/// 1 sign bit, 5 exponent bits, 10 mantissa bits.
fn f16_to_f32(bits: u16) -> f32 {
    let sign = (bits >> 15) & 0x1;
    let exp = (bits >> 10) & 0x1F;
    let mant = bits & 0x3FF;

    let sign_f = if sign == 1 { -1.0 } else { 1.0 };

    if exp == 0 {
        // Subnormal or zero.
        if mant == 0 {
            return sign_f * 0.0;
        }
        let m = mant as f32 / 1024.0;
        return sign_f * m * 2.0_f32.powi(-14);
    }
    if exp == 0x1F {
        // Inf / NaN — clamp to 0 so we never feed garbage into the EMA.
        return 0.0;
    }
    let e = exp as i32 - 15;
    let m = 1.0 + (mant as f32 / 1024.0);
    sign_f * m * 2.0_f32.powi(e)
}

// ── WGSL ─────────────────────────────────────────────────────────────────────

const REDUCER_SHADER: &str = r#"
struct InitUniform {
    luma_weights: vec3<f32>,
    p_norm:       f32,
};

@group(0) @binding(0) var src_tex  : texture_2d<f32>;
@group(0) @binding(1) var src_samp : sampler;
@group(0) @binding(2) var<uniform> u_init: InitUniform;

struct VSOut {
    @builtin(position) pos: vec4<f32>,
    @location(0)       uv:  vec2<f32>,
};

/// Fullscreen triangle — identical to the tonemap's. Three vertices
/// that extend past NDC on two sides; the rasteriser clips the
/// overshoot, leaving a single triangle that covers the whole target.
@vertex
fn vs_fullscreen(@builtin(vertex_index) vi: u32) -> VSOut {
    var pos = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 3.0, -1.0),
        vec2<f32>(-1.0,  3.0),
    );
    var uv = array<vec2<f32>, 3>(
        vec2<f32>(0.0, 1.0),
        vec2<f32>(2.0, 1.0),
        vec2<f32>(0.0, -1.0),
    );
    var out: VSOut;
    out.pos = vec4<f32>(pos[vi], 0.0, 1.0);
    out.uv  = uv[vi];
    return out;
}

/// Init pass: sample HDR, compute L = dot(rgb, weights), emit L^p.
///
/// Clamping `L` to >= 0 before `pow` avoids NaNs when the HDR texture
/// holds a rare negative value (floating-point drift in blends can
/// produce tiny subzero components that leak through the blending
/// math — the single clamp kills them before `pow` goes undefined).
@fragment
fn fs_init(in: VSOut) -> @location(0) vec4<f32> {
    let rgb = textureSample(src_tex, src_samp, in.uv).rgb;
    let l   = max(0.0, dot(rgb, u_init.luma_weights));
    let v   = pow(l, u_init.p_norm);
    return vec4<f32>(v, 0.0, 0.0, 1.0);
}

/// Reduce pass: average the 2×2 neighbourhood of the previous mip.
///
/// Bilinear sampling at the destination texel centre fetches a four-
/// sample average in a single tap — the GPU texture unit does the
/// weighting for us. `src_tex` is bound to the previous mip only;
/// the `reduce_bgl` layout omits the init uniform because this pass
/// needs no per-frame constants.
@fragment
fn fs_reduce(in: VSOut) -> @location(0) vec4<f32> {
    let v = textureSample(src_tex, src_samp, in.uv).r;
    return vec4<f32>(v, 0.0, 0.0, 1.0);
}
"#;

#[cfg(test)]
mod shader_tests {
    #[test]
    fn reducer_shader_validates() {
        crate::render::validate_wgsl("luminance_reducer", super::REDUCER_SHADER);
    }
}
