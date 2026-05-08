//! Bloom pass over the luminous HDR plane.
//!
//! Three render-target fragment passes:
//!
//! 1. **Threshold + cap + downsample** — reads the luminous plane,
//!    discards fragments below the threshold luma, clamps the rest to a
//!    cap (keeps the blur kernel out of `inf` territory when a star is
//!    `1e7+` bright), and writes into a quarter-resolution texture.
//! 2. **Horizontal Gaussian blur** — 9-tap separable kernel.
//! 3. **Vertical Gaussian blur** — same kernel.
//!
//! The composite pass samples the final blurred texture with bilinear
//! filtering so the upscale to full canvas resolution stays smooth.
//!
//! # Format
//!
//! The bloom textures use the same `Rgba16Float` HDR format as the
//! scene targets — the cap (`BLOOM_CAP`) holds blurred values inside the
//! representable range without losing the colour information ACES
//! eventually compresses.

use bytemuck::{Pod, Zeroable};
use std::mem::size_of;

use crate::render::hdr::{HDR_FORMAT, HdrTarget};

const BLOOM_DOWNSCALE: u32 = 4;

/// Luma threshold (Rec. 709) below which a fragment doesn't contribute
/// to the bloom buffer. Pixels with `luma < 1.0` are already in the
/// linear-low region of ACES; their halo would be invisible anyway.
const BLOOM_THRESHOLD: f32 = 1.0;

/// Per-channel clamp on the source pixel before blurring. Decouples
/// halo intensity from raw flux: a star `100×` brighter than another
/// gets `~ln(100)` more halo, not `100×`.
const BLOOM_CAP: f32 = 1.0e3;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ThresholdUniform {
    threshold: f32,
    cap: f32,
    _pad: [f32; 2],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct BlurUniform {
    /// Texel offset in source UV coordinates. `(1/w, 0)` for the
    /// horizontal pass, `(0, 1/h)` for the vertical pass.
    direction: [f32; 2],
    _pad: [f32; 2],
}

pub struct BloomPipeline {
    threshold_pipeline: wgpu::RenderPipeline,
    blur_pipeline: wgpu::RenderPipeline,

    threshold_bgl: wgpu::BindGroupLayout,
    blur_bgl: wgpu::BindGroupLayout,

    threshold_uniform_buf: wgpu::Buffer,
    blur_h_uniform_buf: wgpu::Buffer,
    blur_v_uniform_buf: wgpu::Buffer,

    sampler: wgpu::Sampler,

    /// Ping-pong textures at `BLOOM_DOWNSCALE` smaller than the canvas.
    /// `tex_a` receives the threshold output and the vertical blur
    /// (final); `tex_b` receives the horizontal blur.
    tex_a: Option<wgpu::Texture>,
    view_a: Option<wgpu::TextureView>,
    tex_b: Option<wgpu::Texture>,
    view_b: Option<wgpu::TextureView>,
    current_size: [u32; 2],

    threshold_bg: Option<wgpu::BindGroup>,
    blur_h_bg: Option<wgpu::BindGroup>,
    blur_v_bg: Option<wgpu::BindGroup>,
    bound_hdr_generation: u64,
}

impl BloomPipeline {
    pub fn new(device: &wgpu::Device) -> Self {
        let threshold_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bloom::threshold_bgl"),
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

        let blur_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bloom::blur_bgl"),
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

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("bloom::shader"),
            source: wgpu::ShaderSource::Wgsl(BLOOM_SHADER.into()),
        });

        let threshold_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("bloom::threshold_layout"),
            bind_group_layouts: &[Some(&threshold_bgl)],
            immediate_size: 0,
        });
        let blur_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("bloom::blur_layout"),
            bind_group_layouts: &[Some(&blur_bgl)],
            immediate_size: 0,
        });

        let threshold_pipeline = make_pipeline(
            device,
            &shader,
            &threshold_layout,
            "bloom::threshold_pipeline",
            "fs_threshold",
        );
        let blur_pipeline =
            make_pipeline(device, &shader, &blur_layout, "bloom::blur_pipeline", "fs_blur");

        let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("bloom::sampler"),
            address_mode_u: wgpu::AddressMode::ClampToEdge,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            address_mode_w: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            mipmap_filter: wgpu::MipmapFilterMode::Nearest,
            ..Default::default()
        });

        let threshold_uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("bloom::threshold_uniform"),
            size: size_of::<ThresholdUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let blur_h_uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("bloom::blur_h_uniform"),
            size: size_of::<BlurUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let blur_v_uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("bloom::blur_v_uniform"),
            size: size_of::<BlurUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self {
            threshold_pipeline,
            blur_pipeline,
            threshold_bgl,
            blur_bgl,
            threshold_uniform_buf,
            blur_h_uniform_buf,
            blur_v_uniform_buf,
            sampler,
            tex_a: None,
            view_a: None,
            tex_b: None,
            view_b: None,
            current_size: [0, 0],
            threshold_bg: None,
            blur_h_bg: None,
            blur_v_bg: None,
            bound_hdr_generation: 0,
        }
    }

    fn ensure_textures(&mut self, device: &wgpu::Device, hdr_size: [u32; 2]) -> bool {
        let w = (hdr_size[0] / BLOOM_DOWNSCALE).max(1);
        let h = (hdr_size[1] / BLOOM_DOWNSCALE).max(1);
        if self.tex_a.is_some() && self.current_size == [w, h] {
            return false;
        }
        let make = |label: &str| -> (wgpu::Texture, wgpu::TextureView) {
            let texture = device.create_texture(&wgpu::TextureDescriptor {
                label: Some(label),
                size: wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format: HDR_FORMAT,
                usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                    | wgpu::TextureUsages::TEXTURE_BINDING,
                view_formats: &[],
            });
            let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
            (texture, view)
        };
        let (tex_a, view_a) = make("bloom::tex_a");
        let (tex_b, view_b) = make("bloom::tex_b");
        self.tex_a = Some(tex_a);
        self.view_a = Some(view_a);
        self.tex_b = Some(tex_b);
        self.view_b = Some(view_b);
        self.current_size = [w, h];
        true
    }

    fn rebuild_bind_groups(&mut self, device: &wgpu::Device, hdr: &HdrTarget) {
        let (Some(view_a), Some(view_b)) = (self.view_a.as_ref(), self.view_b.as_ref()) else {
            return;
        };

        self.threshold_bg = Some(device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bloom::threshold_bg"),
            layout: &self.threshold_bgl,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(hdr.view_l()),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::Sampler(&self.sampler),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.threshold_uniform_buf.as_entire_binding(),
                },
            ],
        }));
        self.blur_h_bg = Some(device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bloom::blur_h_bg"),
            layout: &self.blur_bgl,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(view_a),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::Sampler(&self.sampler),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.blur_h_uniform_buf.as_entire_binding(),
                },
            ],
        }));
        self.blur_v_bg = Some(device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bloom::blur_v_bg"),
            layout: &self.blur_bgl,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(view_b),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::Sampler(&self.sampler),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.blur_v_uniform_buf.as_entire_binding(),
                },
            ],
        }));
        self.bound_hdr_generation = hdr.generation();
    }

    pub fn refresh(&mut self, device: &wgpu::Device, hdr: &HdrTarget) {
        let resized = self.ensure_textures(device, hdr.size());
        if resized || self.threshold_bg.is_none() || self.bound_hdr_generation != hdr.generation() {
            self.rebuild_bind_groups(device, hdr);
        }
    }

    pub fn encode(
        &mut self,
        queue: &wgpu::Queue,
        encoder: &mut wgpu::CommandEncoder,
        _hdr: &HdrTarget,
    ) {
        let (Some(view_a), Some(view_b)) = (self.view_a.as_ref(), self.view_b.as_ref()) else {
            return;
        };
        let (Some(threshold_bg), Some(blur_h_bg), Some(blur_v_bg)) =
            (self.threshold_bg.as_ref(), self.blur_h_bg.as_ref(), self.blur_v_bg.as_ref())
        else {
            return;
        };

        queue.write_buffer(
            &self.threshold_uniform_buf,
            0,
            bytemuck::bytes_of(&ThresholdUniform {
                threshold: BLOOM_THRESHOLD,
                cap: BLOOM_CAP,
                _pad: [0.0; 2],
            }),
        );
        let inv_w = 1.0 / self.current_size[0] as f32;
        let inv_h = 1.0 / self.current_size[1] as f32;
        queue.write_buffer(
            &self.blur_h_uniform_buf,
            0,
            bytemuck::bytes_of(&BlurUniform { direction: [inv_w, 0.0], _pad: [0.0; 2] }),
        );
        queue.write_buffer(
            &self.blur_v_uniform_buf,
            0,
            bytemuck::bytes_of(&BlurUniform { direction: [0.0, inv_h], _pad: [0.0; 2] }),
        );

        // Pass 1: HDR_L → tex_a (threshold + cap + downsample by 4×).
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("bloom::threshold_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: view_a,
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
            pass.set_pipeline(&self.threshold_pipeline);
            pass.set_bind_group(0, threshold_bg, &[]);
            pass.draw(0..3, 0..1);
        }

        // Pass 2: tex_a → tex_b (horizontal blur).
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("bloom::blur_h_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: view_b,
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
            pass.set_pipeline(&self.blur_pipeline);
            pass.set_bind_group(0, blur_h_bg, &[]);
            pass.draw(0..3, 0..1);
        }

        // Pass 3: tex_b → tex_a (vertical blur, final).
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("bloom::blur_v_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: view_a,
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
            pass.set_pipeline(&self.blur_pipeline);
            pass.set_bind_group(0, blur_v_bg, &[]);
            pass.draw(0..3, 0..1);
        }
    }

    /// Final blurred bloom view. Sampled by the composite pass.
    pub fn final_view(&self) -> Option<&wgpu::TextureView> {
        self.view_a.as_ref()
    }
}

fn make_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    layout: &wgpu::PipelineLayout,
    label: &str,
    fs_entry: &str,
) -> wgpu::RenderPipeline {
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some(label),
        layout: Some(layout),
        vertex: wgpu::VertexState {
            module: shader,
            entry_point: Some("vs_fullscreen"),
            buffers: &[],
            compilation_options: Default::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: shader,
            entry_point: Some(fs_entry),
            targets: &[Some(wgpu::ColorTargetState {
                format: HDR_FORMAT,
                blend: None,
                write_mask: wgpu::ColorWrites::ALL,
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
    })
}

const BLOOM_SHADER: &str = r#"
struct ThresholdUniform {
    threshold: f32,
    cap:       f32,
    _pad0:     f32,
    _pad1:     f32,
};

struct BlurUniform {
    direction: vec2<f32>,
    _pad0:     f32,
    _pad1:     f32,
};

@group(0) @binding(0) var src_tex  : texture_2d<f32>;
@group(0) @binding(1) var src_samp : sampler;
@group(0) @binding(2) var<uniform> u_threshold: ThresholdUniform;

struct VSOut {
    @builtin(position) pos: vec4<f32>,
    @location(0)       uv:  vec2<f32>,
};

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

const LUMA_WEIGHTS: vec3<f32> = vec3<f32>(0.2126, 0.7152, 0.0722);

@fragment
fn fs_threshold(in: VSOut) -> @location(0) vec4<f32> {
    let rgb = textureSample(src_tex, src_samp, in.uv).rgb;
    let l = dot(rgb, LUMA_WEIGHTS);
    if l < u_threshold.threshold { return vec4<f32>(0.0); }
    let capped = min(rgb, vec3<f32>(u_threshold.cap));
    return vec4<f32>(capped, 1.0);
}

@group(0) @binding(0) var blur_tex  : texture_2d<f32>;
@group(0) @binding(1) var blur_samp : sampler;
@group(0) @binding(2) var<uniform> u_blur: BlurUniform;

/// 9-tap binomial Gaussian, separable. Weights from Pascal's row 8
/// normalised to unit sum: [1, 8, 28, 56, 70, 56, 28, 8, 1] / 256.
@fragment
fn fs_blur(in: VSOut) -> @location(0) vec4<f32> {
    let w0 = 70.0 / 256.0;
    let w1 = 56.0 / 256.0;
    let w2 = 28.0 / 256.0;
    let w3 =  8.0 / 256.0;
    let w4 =  1.0 / 256.0;
    let d  = u_blur.direction;

    var acc = textureSample(blur_tex, blur_samp, in.uv).rgb * w0;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv + d * 1.0).rgb * w1;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv - d * 1.0).rgb * w1;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv + d * 2.0).rgb * w2;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv - d * 2.0).rgb * w2;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv + d * 3.0).rgb * w3;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv - d * 3.0).rgb * w3;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv + d * 4.0).rgb * w4;
    acc = acc + textureSample(blur_tex, blur_samp, in.uv - d * 4.0).rgb * w4;
    return vec4<f32>(acc, 1.0);
}
"#;
