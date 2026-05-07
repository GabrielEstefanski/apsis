//! Composite + tonemap pass: dual HDR → sRGB swapchain.
//!
//! Samples the reflective and luminous planes from
//! [`crate::render::hdr::HdrTarget`], multiplies the reflective sample by
//! the exposure scalar (auto-exposure × user EV), sums both in HDR, and
//! runs a single ACES filmic curve.
//!
//! # Gamma
//!
//! No manual gamma correction in this shader. The swapchain is assumed to be
//! an sRGB format, in which case the GPU applies the OETF on write. If the
//! swapchain is a non-sRGB format, colours come out looking darker — in that
//! case the caller should either request an sRGB surface or encode gamma
//! explicitly at the tonemap output.
//!
//! # Resize tracking
//!
//! The bind group samples both HDR texture views, which are recreated on
//! every canvas resize. [`TonemapPipeline::refresh_if_resized`] checks the
//! HDR target's generation counter and rebuilds the bind group lazily.

use std::mem::size_of;

use bytemuck::{Pod, Zeroable};

use crate::render::hdr::HdrTarget;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct TonemapUniform {
    /// Linear multiplier on the reflective plane.
    exposure_r: f32,
    _pad: [f32; 3],
}

impl Default for TonemapUniform {
    fn default() -> Self {
        Self { exposure_r: 1.0, _pad: [0.0; 3] }
    }
}

pub struct TonemapPipeline {
    pipeline: wgpu::RenderPipeline,
    sampler: wgpu::Sampler,
    uniform_buf: wgpu::Buffer,
    bgl: wgpu::BindGroupLayout,
    bind_group: Option<wgpu::BindGroup>,
    /// Last observed [`HdrTarget::generation`] for the bind group — when it
    /// disagrees with the target's current generation, the bind group is
    /// rebuilt.
    bound_generation: u64,
    exposure: f32,
}

impl TonemapPipeline {
    pub fn new(device: &wgpu::Device, swapchain_format: wgpu::TextureFormat) -> Self {
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("tonemap::bgl"),
            entries: &[
                // 0: reflective HDR
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
                // 1: luminous HDR
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // 2: bloom (blurred luminous, bilinear-upscaled at sample time)
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // 3: shared sampler
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
                // 4: tonemap uniform
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
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
            label: Some("tonemap::shader"),
            source: wgpu::ShaderSource::Wgsl(TONEMAP_SHADER.into()),
        });

        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("tonemap::layout"),
            bind_group_layouts: &[Some(&bgl)],
            immediate_size: 0,
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("tonemap::pipeline"),
            layout: Some(&layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_fullscreen"),
                buffers: &[],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_tonemap"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: swapchain_format,
                    // Composite the tonemapped scene on top of whatever egui
                    // drew behind the callback rect. The HDR target was
                    // cleared transparent, so outside the scene its alpha is
                    // zero and the backdrop shows through.
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
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
        });

        let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("tonemap::sampler"),
            address_mode_u: wgpu::AddressMode::ClampToEdge,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            address_mode_w: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            mipmap_filter: wgpu::MipmapFilterMode::Nearest,
            ..Default::default()
        });

        let uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("tonemap::uniform"),
            size: size_of::<TonemapUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self {
            pipeline,
            sampler,
            uniform_buf,
            bgl,
            bind_group: None,
            bound_generation: 0,
            exposure: 1.0,
        }
    }

    pub fn set_exposure(&mut self, exposure: f32) {
        self.exposure = exposure.max(0.0);
    }

    /// Rebuilds the sampled-texture bind group when the HDR target has been
    /// reallocated since the last call.
    pub fn refresh_if_resized(
        &mut self,
        device: &wgpu::Device,
        hdr: &HdrTarget,
        bloom_view: &wgpu::TextureView,
    ) {
        if self.bind_group.is_some() && self.bound_generation == hdr.generation() {
            return;
        }
        self.bind_group = Some(device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("tonemap::bg"),
            layout: &self.bgl,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(hdr.view_r()),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(hdr.view_l()),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(bloom_view),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: wgpu::BindingResource::Sampler(&self.sampler),
                },
                wgpu::BindGroupEntry { binding: 4, resource: self.uniform_buf.as_entire_binding() },
            ],
        }));
        self.bound_generation = hdr.generation();
    }

    pub fn upload(&self, queue: &wgpu::Queue) {
        let u = TonemapUniform { exposure_r: self.exposure, _pad: [0.0; 3] };
        queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&u));
    }

    pub fn draw(&self, pass: &mut wgpu::RenderPass<'_>) {
        let Some(bg) = self.bind_group.as_ref() else {
            return;
        };
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, bg, &[]);
        pass.draw(0..3, 0..1);
    }
}

// ── WGSL ──────────────────────────────────────────────────────────────────────

const TONEMAP_SHADER: &str = r#"
struct TonemapUniform {
    exposure_r: f32,
    _pad0:      f32,
    _pad1:      f32,
    _pad2:      f32,
};

@group(0) @binding(0) var hdr_r_tex : texture_2d<f32>;
@group(0) @binding(1) var hdr_l_tex : texture_2d<f32>;
@group(0) @binding(2) var bloom_tex : texture_2d<f32>;
@group(0) @binding(3) var hdr_samp  : sampler;
@group(0) @binding(4) var<uniform> u: TonemapUniform;

struct VSOut {
    @builtin(position) pos: vec4<f32>,
    @location(0)       uv:  vec2<f32>,
};

/// Fullscreen triangle — three vertices cover the whole viewport with no
/// vertex buffer. Simpler than a quad and avoids a diagonal seam in the UVs.
@vertex
fn vs_fullscreen(@builtin(vertex_index) vi: u32) -> VSOut {
    // Three points that span NDC [-1, +1] × [-1, +1] with the unused corner
    // clipped off. UVs matched so sampling lands on [0, 1] within the viewport.
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

/// ACES filmic tonemap (Krzysztof Narkowicz's fit).
///
/// Closed-form rational approximation of the full ACES RRT + ODT pipeline.
/// Fast, stable across the HDR range, and preserves highlight hue better
/// than Reinhard when multiple bright sources share a pixel.
fn aces(x: vec3<f32>) -> vec3<f32> {
    let a = 2.51;
    let b = 0.03;
    let c = 2.43;
    let d = 0.59;
    let e = 0.14;
    return clamp((x * (a * x + b)) / (x * (c * x + d) + e),
                 vec3<f32>(0.0),
                 vec3<f32>(1.0));
}

@fragment
fn fs_tonemap(in: VSOut) -> @location(0) vec4<f32> {
    let src_r = textureSample(hdr_r_tex, hdr_samp, in.uv);
    let src_l = textureSample(hdr_l_tex, hdr_samp, in.uv);
    let bloom = textureSample(bloom_tex, hdr_samp, in.uv);
    let combined = src_r.rgb * u.exposure_r + src_l.rgb + bloom.rgb;
    let lit = aces(combined);
    let alpha = max(max(src_r.a, src_l.a), bloom.a);
    return vec4<f32>(lit, alpha);
}
"#;

#[cfg(test)]
mod shader_tests {
    #[test]
    fn tonemap_shader_validates() {
        crate::render::validate_wgsl("tonemap", super::TONEMAP_SHADER);
    }
}
