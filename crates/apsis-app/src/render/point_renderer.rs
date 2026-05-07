//! Point-source rendering for sub-pixel reflective bodies.
//!
//! When a body's projected radius drops below ~1 px, the disc impostor
//! either fails to rasterise or aliases hard. Each such body is drawn
//! instead as a single Gaussian sprite — 5×5 px support, σ ≈ 0.8 px,
//! normalised so `Σ kernel = 1`. The `intensity_linear` value passed in
//! is the body's full HDR contribution; the kernel just spreads it
//! across the support window.
//!
//! Pre-multiplied additive blend. Targets the reflective HDR plane;
//! luminous bodies use the disc path with a min-pixel floor and feed
//! the bloom pass instead.

use bytemuck::{Pod, Zeroable};
use std::mem::size_of;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
pub struct PointInstance {
    pub screen_pos: [f32; 2],
    pub intensity: f32,
    pub _pad0: f32,
    pub color: [f32; 3],
    pub _pad1: f32,
}

impl PointInstance {
    pub fn new(screen_pos: [f32; 2], intensity: f32, color: [f32; 3]) -> Self {
        Self { screen_pos, intensity, _pad0: 0.0, color, _pad1: 0.0 }
    }
}

pub struct PointRenderer {
    pipeline: wgpu::RenderPipeline,
    instance_buf: wgpu::Buffer,
    instance_cap: u32,
}

impl PointRenderer {
    pub fn new(
        device: &wgpu::Device,
        screen_bgl: &wgpu::BindGroupLayout,
        format: wgpu::TextureFormat,
    ) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("point::shader"),
            source: wgpu::ShaderSource::Wgsl(POINT_SHADER.into()),
        });

        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("point::layout"),
            bind_group_layouts: &[Some(screen_bgl)],
            immediate_size: 0,
        });

        let attrs = wgpu::vertex_attr_array![
            0 => Float32x2,
            1 => Float32,
            2 => Float32,
            3 => Float32x3,
            4 => Float32,
        ];

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("point::pipeline"),
            layout: Some(&layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_point"),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: size_of::<PointInstance>() as u64,
                    step_mode: wgpu::VertexStepMode::Instance,
                    attributes: &attrs,
                }],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_point"),
                targets: &[Some(wgpu::ColorTargetState {
                    format,
                    blend: Some(wgpu::BlendState {
                        color: wgpu::BlendComponent {
                            src_factor: wgpu::BlendFactor::One,
                            dst_factor: wgpu::BlendFactor::One,
                            operation: wgpu::BlendOperation::Add,
                        },
                        alpha: wgpu::BlendComponent {
                            src_factor: wgpu::BlendFactor::One,
                            dst_factor: wgpu::BlendFactor::One,
                            operation: wgpu::BlendOperation::Add,
                        },
                    }),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState::default(),
            depth_stencil: None,
            multisample: Default::default(),
            multiview_mask: None,
            cache: None,
        });

        let instance_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("point::instances"),
            size: (size_of::<PointInstance>() as u64).max(64),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self { pipeline, instance_buf, instance_cap: 1 }
    }

    pub fn upload(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        instances: &[PointInstance],
    ) -> u32 {
        let n = instances.len() as u32;
        if n == 0 {
            return 0;
        }
        if n > self.instance_cap {
            let cap = (n * 2).max(256);
            self.instance_buf = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("point::instances"),
                size: cap as u64 * size_of::<PointInstance>() as u64,
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.instance_cap = cap;
        }
        queue.write_buffer(&self.instance_buf, 0, bytemuck::cast_slice(instances));
        n
    }

    pub fn draw(
        &self,
        pass: &mut wgpu::RenderPass<'_>,
        screen_bg: &wgpu::BindGroup,
        instance_count: u32,
    ) {
        if instance_count == 0 {
            return;
        }
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, screen_bg, &[]);
        pass.set_vertex_buffer(0, self.instance_buf.slice(..));
        pass.draw(0..6, 0..instance_count);
    }
}

const POINT_SHADER: &str = r#"
struct ScreenUniform {
    size:         vec2<f32>,
    viewport_min: vec2<f32>,
};

@group(0) @binding(0) var<uniform> screen: ScreenUniform;

struct PointInstance {
    @location(0) screen_pos: vec2<f32>,
    @location(1) intensity:  f32,
    @location(2) _pad0:      f32,
    @location(3) color:      vec3<f32>,
    @location(4) _pad1:      f32,
};

struct VSOut {
    @builtin(position) pos:       vec4<f32>,
    @location(0)       color:     vec3<f32>,
    @location(1)       intensity: f32,
    @location(2)       offset_px: vec2<f32>,
};

const HALF_EXTENT_PX: f32 = 2.5;
const KERNEL_SIGMA_PX: f32 = 0.8;
const KERNEL_NORM: f32 = 1.0 / (2.0 * 3.141592653589793 * KERNEL_SIGMA_PX * KERNEL_SIGMA_PX);

fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let local = p - screen.viewport_min;
    let x =  (local.x / screen.size.x) * 2.0 - 1.0;
    let y = -(local.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

@vertex
fn vs_point(
    @builtin(vertex_index) vi: u32,
    instance: PointInstance,
) -> VSOut {
    var corners = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );
    let corner = corners[vi];
    let offset_px = corner * HALF_EXTENT_PX;
    let pixel_pos = instance.screen_pos + offset_px;

    var out: VSOut;
    out.pos       = to_ndc(pixel_pos);
    out.color     = instance.color;
    out.intensity = instance.intensity;
    out.offset_px = offset_px;
    return out;
}

@fragment
fn fs_point(in: VSOut) -> @location(0) vec4<f32> {
    let r2 = dot(in.offset_px, in.offset_px);
    let kernel = exp(-r2 / (2.0 * KERNEL_SIGMA_PX * KERNEL_SIGMA_PX)) * KERNEL_NORM;
    let amount = in.intensity * kernel;
    if amount <= 0.0 { discard; }
    return vec4<f32>(in.color * amount, amount);
}
"#;
