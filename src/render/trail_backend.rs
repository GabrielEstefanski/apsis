//! GPU trail renderer using a persistent ring-buffer on the GPU.
//!
//! # Overview
//!
//! [`TrailRenderer::submit`] is called once per frame while trails are visible.
//! It drains the [`TrailBuffer`]'s dirty flags, snapshots only the changed data,
//! and registers a [`TrailCallback`] with egui's WGPU integration.
//!
//! The callback's `prepare` phase uploads the minimal set of bytes to the GPU
//! (one column = N × 8 bytes on the common path) and updates the per-frame
//! state uniform.  The `paint` phase issues a single `draw` call:
//!
//! ```text
//! draw(vertices = 6 × (cap − 1), instances = n_bodies)
//! ```
//!
//! # Buffer layout
//!
//! | Buffer | Binding | Usage |
//! |--------|---------|-------|
//! | `screen_buf` | group 0, binding 0 | `ScreenUniform` — logical viewport size |
//! | `pos_buf`    | group 1, binding 0 | `array<vec2<f32>>` — column-major positions |
//! | `color_buf`  | group 1, binding 1 | `array<vec4<f32>>` — per-body RGBA |
//! | `state_buf`  | group 1, binding 2 | `GpuTrailState` — ring-buffer parameters |
//!
//! # Vertex shader
//!
//! Each vertex is identified by `(vertex_index, instance_index)`:
//! - `instance_index` = body index
//! - `seg_idx = vertex_index / 6` — segment along the trail (0 = newest)
//! - `corner   = vertex_index % 6` — corner of the screen-aligned quad
//!
//! Segments older than `trail_len − 1` or containing NaN positions produce
//! degenerate triangles (clip-space position `(2,2,0,1)`) that are invisible.

use std::mem::size_of;
use std::sync::atomic::{AtomicU32, Ordering};

use bytemuck::{Pod, Zeroable};
use eframe::egui;
use eframe::egui_wgpu::{self, CallbackTrait, ScreenDescriptor};

use crate::core::trail_buffer::{PendingPositions, TrailBuffer};

// ── GPU-mapped data types ─────────────────────────────────────────────────────

/// Screen dimensions uniform (matches the layout in `wgpu_backend.rs`).
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ScreenUniform {
    size: [f32; 2],
    _pad: [f32; 2],
}

impl ScreenUniform {
    fn from_descriptor(sd: &ScreenDescriptor) -> Self {
        Self {
            size: [
                sd.size_in_pixels[0] as f32 / sd.pixels_per_point,
                sd.size_in_pixels[1] as f32 / sd.pixels_per_point,
            ],
            _pad: [0.0; 2],
        }
    }
}

/// Per-frame trail state sent to the vertex shader.
///
/// Layout is exactly 32 bytes, meeting the `minUniformBufferOffsetAlignment`
/// requirement of all common WGPU back-ends.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct GpuTrailState {
    /// Ring-buffer write head (= index of the *next* column to be written).
    head: u32,
    /// Number of valid columns (capped at `cap`).
    trail_len: u32,
    /// Body count — matches the inner dimension of the position matrix.
    n_bodies: u32,
    /// Ring-buffer depth (outer dimension of the position matrix).
    cap: u32,
    /// Canvas centre in logical screen pixels.
    center: [f32; 2],
    /// Pixels per world unit.
    scale: f32,
    _pad: f32,
}

// ── Persistent GPU resources ──────────────────────────────────────────────────

/// GPU-side state that persists across frames.
///
/// Stored in [`egui_wgpu::CallbackResources`] and initialised lazily.
/// `pos_buf` and `color_buf` are reallocated whenever the body count or
/// ring-buffer depth changes; all other resources remain stable.
struct GpuTrailResources {
    pipeline: wgpu::RenderPipeline,
    screen_buf: wgpu::Buffer,
    state_buf: wgpu::Buffer,
    pos_buf: wgpu::Buffer,
    color_buf: wgpu::Buffer,
    bind_group_screen: wgpu::BindGroup,
    bind_group_trail: wgpu::BindGroup,
    /// Kept alive for realloc (recreating `bind_group_trail`).
    trail_bgl: wgpu::BindGroupLayout,
    /// Cached dimensions; used to detect topology changes.
    n_bodies: u32,
    cap: u32,
}

impl GpuTrailResources {
    /// Creates all GPU resources, using 1-element dummy buffers for the
    /// variable-size `pos_buf` / `color_buf`.  [`realloc`](Self::realloc)
    /// is called immediately on the first frame with real dimensions.
    fn init(device: &wgpu::Device) -> Self {
        let screen_bgl = build_screen_bgl(device);
        let trail_bgl = build_trail_bgl(device);
        let shader = build_trail_shader(device);
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("trail::pipeline_layout"),
            bind_group_layouts: &[Some(&screen_bgl), Some(&trail_bgl)],
            immediate_size: 0,
        });
        let pipeline = build_trail_pipeline(device, &shader, &pipeline_layout);

        let screen_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::screen"),
            size: size_of::<ScreenUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let state_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::state"),
            size: size_of::<GpuTrailState>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group_screen = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("trail::bg_screen"),
            layout: &screen_bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: screen_buf.as_entire_binding(),
            }],
        });

        // Dummy 1-element buffers — replaced on first realloc.
        let pos_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::pos_dummy"),
            size: 8,  // one vec2<f32>
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let color_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::color_dummy"),
            size: 16,  // one vec4<f32>
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group_trail = build_trail_bind_group(
            device, &trail_bgl, &pos_buf, &color_buf, &state_buf,
        );

        Self {
            pipeline,
            screen_buf,
            state_buf,
            pos_buf,
            color_buf,
            bind_group_screen,
            bind_group_trail,
            trail_bgl,
            n_bodies: 0,
            cap: 0,
        }
    }

    /// Reallocates `pos_buf` and `color_buf` for a new `(n_bodies, cap)`.
    ///
    /// Recreates `bind_group_trail` to point at the new buffers.
    /// All other resources are unaffected.
    fn realloc(&mut self, device: &wgpu::Device, n_bodies: u32, cap: u32) {
        let pos_bytes = (n_bodies as u64 * cap as u64 * 8).max(8);
        let color_bytes = (n_bodies as u64 * 16).max(16);

        self.pos_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::positions"),
            size: pos_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.color_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::colors"),
            size: color_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.bind_group_trail = build_trail_bind_group(
            device,
            &self.trail_bgl,
            &self.pos_buf,
            &self.color_buf,
            &self.state_buf,
        );
        self.n_bodies = n_bodies;
        self.cap = cap;
    }
}

// ── Per-frame callback ────────────────────────────────────────────────────────

/// Per-frame trail paint callback registered with egui's WGPU integration.
///
/// Carries CPU-side snapshots of changed GPU data.  Persistent GPU resources
/// live in [`GpuTrailResources`] inside [`egui_wgpu::CallbackResources`].
struct TrailCallback {
    state: GpuTrailState,
    n_bodies: u32,
    cap: u32,
    /// Which parts of the position buffer need uploading.
    pending_pos: PendingPositions,
    /// Full position matrix snapshot (Some when `pending_pos` is `Full`).
    full_snapshot: Option<Vec<[f32; 2]>>,
    /// Per-column position snapshots (non-empty when `pending_pos` is `Columns`).
    col_snapshots: Vec<(u32, Vec<[f32; 2]>)>,
    /// Colour data (Some when colours are dirty).
    colors: Option<Vec<[f32; 4]>>,
    /// Written by `prepare`, read by `paint`.
    vertex_count: AtomicU32,
    /// Written by `prepare`, read by `paint`.
    instance_count: AtomicU32,
}

impl CallbackTrait for TrailCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        screen_descriptor: &ScreenDescriptor,
        _encoder: &mut wgpu::CommandEncoder,
        resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        // Lazy init.
        if resources.get::<GpuTrailResources>().is_none() {
            resources.insert(GpuTrailResources::init(device));
        }
        let gpu = resources
            .get_mut::<GpuTrailResources>()
            .expect("GpuTrailResources present after lazy init");

        // Realloc when the body count or ring-buffer depth changes.
        if gpu.n_bodies != self.n_bodies || gpu.cap != self.cap {
            gpu.realloc(device, self.n_bodies, self.cap);
        }

        // Screen uniform — uploaded every frame (16 bytes, negligible cost).
        queue.write_buffer(
            &gpu.screen_buf,
            0,
            bytemuck::bytes_of(&ScreenUniform::from_descriptor(screen_descriptor)),
        );

        // State uniform — always changes (head, trail_len, center, scale).
        queue.write_buffer(&gpu.state_buf, 0, bytemuck::bytes_of(&self.state));

        // Position upload — incremental on the common path.
        match &self.pending_pos {
            PendingPositions::Full => {
                if let Some(data) = &self.full_snapshot {
                    queue.write_buffer(&gpu.pos_buf, 0, bytemuck::cast_slice(data));
                }
            }
            PendingPositions::Columns(_) => {
                let n = self.n_bodies as u64;
                for (col, data) in &self.col_snapshots {
                    let byte_offset = *col as u64 * n * 8; // 8 bytes per vec2<f32>
                    queue.write_buffer(&gpu.pos_buf, byte_offset, bytemuck::cast_slice(data));
                }
            }
            PendingPositions::Clean => {}
        }

        // Colour upload — only when a body colour changed.
        if let Some(colors) = &self.colors {
            queue.write_buffer(&gpu.color_buf, 0, bytemuck::cast_slice(colors));
        }

        // Compute draw counts for the paint phase.
        let vtx = if self.state.trail_len >= 2 {
            6 * (self.cap - 1)
        } else {
            0
        };
        self.vertex_count.store(vtx, Ordering::Relaxed);
        self.instance_count.store(self.n_bodies, Ordering::Relaxed);

        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        pass: &mut wgpu::RenderPass<'static>,
        resources: &egui_wgpu::CallbackResources,
    ) {
        let Some(gpu) = resources.get::<GpuTrailResources>() else {
            return;
        };

        // SAFETY: `GpuTrailResources` is stored in `CallbackResources`, which
        // enforces `'static` on all inserted values and is application-scoped.
        // We need `'static` so that the render pass (which borrows the bind
        // groups) can satisfy its own lifetime bound.
        let gpu: &'static GpuTrailResources =
            unsafe { &*(gpu as *const GpuTrailResources) };

        let vtx = self.vertex_count.load(Ordering::Relaxed);
        let inst = self.instance_count.load(Ordering::Relaxed);
        if vtx == 0 || inst == 0 {
            return;
        }

        pass.set_pipeline(&gpu.pipeline);
        pass.set_bind_group(0, &gpu.bind_group_screen, &[]);
        pass.set_bind_group(1, &gpu.bind_group_trail, &[]);
        pass.draw(0..vtx, 0..inst);
    }
}

// ── Public interface ──────────────────────────────────────────────────────────

/// Submits a GPU trail render callback for the current frame.
///
/// Call once per frame while trails are enabled, **before** the body
/// primitive callback so that trails render underneath bodies.
pub struct TrailRenderer;

impl TrailRenderer {
    /// Drains the trail buffer's dirty flags, snapshots any changed data, and
    /// registers a [`TrailCallback`] with egui's WGPU painter.
    ///
    /// # Parameters
    ///
    /// - `ui`: The egui UI for the canvas panel (used to access the painter).
    /// - `rect`: The canvas rectangle (determines the scissor/viewport region).
    /// - `trail_buf`: Mutable reference — dirty flags are drained each call.
    /// - `center`: Canvas centre in logical screen pixels (`[cx, cy]`).
    /// - `scale`: World-to-screen scale factor (logical pixels per world unit).
    pub fn submit(
        ui: &egui::Ui,
        rect: egui::Rect,
        trail_buf: &mut TrailBuffer,
        center: [f32; 2],
        scale: f32,
    ) {
        // Drain flags even if not renderable so they don't accumulate.
        let pending_pos = trail_buf.take_positions_upload();
        let colors_dirty = trail_buf.take_colors_dirty();

        if !trail_buf.is_renderable() {
            return;
        }

        let n_bodies = trail_buf.n_bodies();
        let cap = trail_buf.capacity();

        // Snapshot only the data that changed.
        let full_snapshot = if matches!(pending_pos, PendingPositions::Full) {
            Some(trail_buf.positions().to_vec())
        } else {
            None
        };

        let col_snapshots = if let PendingPositions::Columns(ref cols) = pending_pos {
            cols.iter()
                .map(|&col| (col, trail_buf.column_slice(col).to_vec()))
                .collect()
        } else {
            Vec::new()
        };

        let colors = if colors_dirty {
            Some(trail_buf.colors().to_vec())
        } else {
            None
        };

        let state = GpuTrailState {
            head: trail_buf.head(),
            trail_len: trail_buf.len(),
            n_bodies,
            cap,
            center,
            scale,
            _pad: 0.0,
        };

        ui.painter().add(egui_wgpu::Callback::new_paint_callback(
            rect,
            TrailCallback {
                state,
                n_bodies,
                cap,
                pending_pos,
                full_snapshot,
                col_snapshots,
                colors,
                vertex_count: AtomicU32::new(0),
                instance_count: AtomicU32::new(0),
            },
        ));
    }
}

// ── Pipeline / shader builders ────────────────────────────────────────────────

fn build_screen_bgl(device: &wgpu::Device) -> wgpu::BindGroupLayout {
    device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("trail::screen_bgl"),
        entries: &[wgpu::BindGroupLayoutEntry {
            binding: 0,
            visibility: wgpu::ShaderStages::VERTEX,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        }],
    })
}

fn build_trail_bgl(device: &wgpu::Device) -> wgpu::BindGroupLayout {
    device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("trail::trail_bgl"),
        entries: &[
            // binding 0: position storage buffer
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: true },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            // binding 1: colour storage buffer
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: true },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            // binding 2: trail state uniform
            wgpu::BindGroupLayoutEntry {
                binding: 2,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    })
}

fn build_trail_bind_group(
    device: &wgpu::Device,
    layout: &wgpu::BindGroupLayout,
    pos_buf: &wgpu::Buffer,
    color_buf: &wgpu::Buffer,
    state_buf: &wgpu::Buffer,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("trail::bg_trail"),
        layout,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: pos_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: color_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: state_buf.as_entire_binding() },
        ],
    })
}

fn build_trail_shader(device: &wgpu::Device) -> wgpu::ShaderModule {
    device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("trail::shader"),
        source: wgpu::ShaderSource::Wgsl(TRAIL_SHADER.into()),
    })
}

fn build_trail_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    layout: &wgpu::PipelineLayout,
) -> wgpu::RenderPipeline {
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("trail::pipeline"),
        layout: Some(layout),
        vertex: wgpu::VertexState {
            module: shader,
            entry_point: Some("vs_trail"),
            // All vertex data comes from storage buffers; no vertex buffers needed.
            buffers: &[],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: shader,
            entry_point: Some("fs_trail"),
            targets: &[Some(wgpu::ColorTargetState {
                format: wgpu::TextureFormat::Bgra8Unorm,
                blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                write_mask: wgpu::ColorWrites::ALL,
            })],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::TriangleList,
            ..Default::default()
        },
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        multiview_mask: None,
        cache: None,
    })
}

// ── WGSL shader ───────────────────────────────────────────────────────────────

/// Trail WGSL shader.
///
/// ## Vertex shader
///
/// Procedural — no vertex buffers.  Each draw call covers all bodies
/// (`instance_index = body_idx`) and all segments (`vertex_index / 6 = seg_idx`).
///
/// For each segment the shader looks up two ring-buffer columns, converts their
/// world-space positions to screen-pixel coordinates, and extrudes a
/// screen-aligned quad (6 vertices, 2 triangles) along the segment perpendicular.
///
/// Fade: `t = seg_idx / (trail_len − 1)`, `fade = (1 − t)^2.2`.
/// Newer segments are brighter and slightly wider.
///
/// Invalid segments (NaN position, zero-length direction, or beyond `trail_len`)
/// produce a degenerate triangle at `(2, 2, 0, 1)` which is clipped silently.
const TRAIL_SHADER: &str = r#"
struct ScreenUniform {
    size: vec2<f32>,
    _pad: vec2<f32>,
};

struct TrailState {
    head:      u32,
    trail_len: u32,
    n_bodies:  u32,
    cap:       u32,
    center:    vec2<f32>,
    scale:     f32,
    _pad:      f32,
};

@group(0) @binding(0) var<uniform> screen: ScreenUniform;
@group(1) @binding(0) var<storage, read> positions: array<vec2<f32>>;
@group(1) @binding(1) var<storage, read> colors:    array<vec4<f32>>;
@group(1) @binding(2) var<uniform> state: TrailState;

struct VertOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0)       color:    vec4<f32>,
};

// Maps a logical-pixel position to WGPU clip space (NDC, Y-up).
fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let x =  (p.x / screen.size.x) * 2.0 - 1.0;
    let y = -(p.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

// Converts a world-space position to logical screen pixels.
fn world_to_screen(w: vec2<f32>) -> vec2<f32> {
    return state.center + w * state.scale;
}

// Returns a clipped-away (invisible) vertex.
fn degenerate() -> VertOut {
    var out: VertOut;
    out.clip_pos = vec4<f32>(2.0, 2.0, 0.0, 1.0);
    out.color    = vec4<f32>(0.0);
    return out;
}

// Returns the screen-pixel position of a quad corner.
//
// Corners 0–5 form two CCW triangles (TriangleList):
//   Triangle 1: 0(old,+), 1(new,+), 2(new,-)
//   Triangle 2: 3(old,+), 4(new,-), 5(old,-)
fn corner_pos(
    c:       u32,
    old_pos: vec2<f32>,
    new_pos: vec2<f32>,
    perp:    vec2<f32>,
    hw:      f32,
) -> vec2<f32> {
    if c == 0u || c == 3u {
        return old_pos + perp * hw;
    } else if c == 1u {
        return new_pos + perp * hw;
    } else if c == 2u || c == 4u {
        return new_pos - perp * hw;
    } else {
        return old_pos - perp * hw;
    }
}

@vertex
fn vs_trail(
    @builtin(vertex_index)   vi:       u32,
    @builtin(instance_index) body_idx: u32,
) -> VertOut {
    let seg_idx = vi / 6u;
    let corner  = vi % 6u;

    // Clamp early: only trail_len − 1 segments are renderable.
    if state.trail_len < 2u || seg_idx >= state.trail_len - 1u {
        return degenerate();
    }

    // Ring-buffer column lookup: seg_idx=0 → newest segment.
    let col_new = (state.head + state.cap - 1u - seg_idx) % state.cap;
    let col_old = (state.head + state.cap - 2u - seg_idx) % state.cap;

    let p_new = positions[col_new * state.n_bodies + body_idx];
    let p_old = positions[col_old * state.n_bodies + body_idx];

    // NaN sentinel — unwritten slots.
    if p_new.x != p_new.x || p_old.x != p_old.x {
        return degenerate();
    }

    let s_new = world_to_screen(p_new);
    let s_old = world_to_screen(p_old);

    let delta = s_new - s_old;
    let seg_len = length(delta);
    if seg_len < 1e-4 {
        return degenerate();
    }

    let dir  = delta / seg_len;
    let perp = vec2<f32>(-dir.y, dir.x);

    // Fade from newest (t=0, bright) to oldest (t=1, transparent).
    let t    = f32(seg_idx) / f32(state.trail_len - 1u);
    let fade = pow(1.0 - t, 2.2);
    let hw   = 0.2 + fade * 0.4; // half-width in logical pixels

    let pos = corner_pos(corner, s_old, s_new, perp, hw);
    let col = colors[body_idx];

    var out: VertOut;
    out.clip_pos = to_ndc(pos);
    out.color    = vec4<f32>(col.rgb * fade, col.a * fade);
    return out;
}

@fragment
fn fs_trail(in: VertOut) -> @location(0) vec4<f32> {
    return in.color;
}
"#;
