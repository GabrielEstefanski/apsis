use bytemuck::{Pod, Zeroable};
use std::mem::size_of;

use crate::app::render_hints::BodyRenderHints;
use crate::core::trail_buffer::TrailBuffer;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct TrailState {
    pub head: u32,
    pub trail_len: u32,
    pub n_bodies: u32,
    pub cap: u32,
    pub center: [f32; 2],
    pub scale: f32,
    pub trail_width: f32,
}

pub struct TrailRenderer {
    pipeline: wgpu::RenderPipeline,

    pos_buf: wgpu::Buffer,
    color_buf: wgpu::Buffer,
    state_buf: wgpu::Buffer,

    bind_group: wgpu::BindGroup,
    bind_group_layout: wgpu::BindGroupLayout,

    n_bodies: u32,
    cap: u32,
}

impl TrailRenderer {
    pub fn new(
        device: &wgpu::Device,
        screen_bgl: &wgpu::BindGroupLayout,
        format: wgpu::TextureFormat,
    ) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("trail::bgl"),
            entries: &[
                // positions
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
                // colors
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
                // state
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
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("trail::shader"),
            source: wgpu::ShaderSource::Wgsl(TRAIL_SHADER.into()),
        });

        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("trail::layout"),
            bind_group_layouts: &[Some(screen_bgl), Some(&bind_group_layout)],
            immediate_size: 0,
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("trail::pipeline"),
            layout: Some(&layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_trail"),
                buffers: &[],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_trail"),
                targets: &[Some(wgpu::ColorTargetState {
                    format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
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

        let pos_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::pos"),
            size: 8,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let color_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::color"),
            size: 16,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let state_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trail::state"),
            size: size_of::<TrailState>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("trail::bg"),
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: pos_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: color_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: state_buf.as_entire_binding(),
                },
            ],
        });

        Self {
            pipeline,
            pos_buf,
            color_buf,
            state_buf,
            bind_group,
            bind_group_layout,
            n_bodies: 0,
            cap: 0,
        }
    }

    pub fn upload(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        trail: &TrailBuffer,
        center: [f32; 2],
        scale: f32,
        trail_width: f32,
    ) {
        let n_bodies = trail.n_bodies();
        let cap = trail.capacity();

        if self.n_bodies != n_bodies || self.cap != cap {
            self.n_bodies = n_bodies;
            self.cap = cap;

            let pos_size = (n_bodies * cap * 8).max(8);
            let color_size = (n_bodies * 16).max(16);

            self.pos_buf = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("trail::pos"),
                size: pos_size as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });

            self.color_buf = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("trail::color"),
                size: color_size as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });

            self.bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("trail::bg"),
                layout: &self.bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: self.pos_buf.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: self.color_buf.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: self.state_buf.as_entire_binding(),
                    },
                ],
            });
        }

        queue.write_buffer(&self.pos_buf, 0, bytemuck::cast_slice(trail.positions()));
        queue.write_buffer(&self.color_buf, 0, bytemuck::cast_slice(trail.colors()));

        let state = TrailState {
            head: trail.head(),
            trail_len: trail.len(),
            n_bodies,
            cap,
            center,
            scale,
            trail_width,
        };

        queue.write_buffer(&self.state_buf, 0, bytemuck::bytes_of(&state));
    }

    pub fn draw(&self, pass: &mut wgpu::RenderPass<'_>, screen_bg: &wgpu::BindGroup) {
        if self.cap < 2 || self.n_bodies == 0 {
            return;
        }

        let vtx = 6 * (self.cap - 1);

        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, screen_bg, &[]);
        pass.set_bind_group(1, &self.bind_group, &[]);
        pass.draw(0..vtx, 0..self.n_bodies);
    }

    pub fn bind_group_layout(&self) -> &wgpu::BindGroupLayout {
        &self.bind_group_layout
    }
}

const TRAIL_SHADER: &str = r#"
struct ScreenUniform {
    size: vec2<f32>,
    _pad: vec2<f32>,
};

@group(0) @binding(0)
var<uniform> screen: ScreenUniform;

struct TrailState {
    head: u32,
    trail_len: u32,
    n_bodies: u32,
    cap: u32,
    center: vec2<f32>,
    scale: f32,
    trail_width: f32,
};

@group(1) @binding(0)
var<storage, read> positions: array<vec2<f32>>;

@group(1) @binding(1)
var<storage, read> colors: array<vec4<f32>>;

@group(1) @binding(2)
var<uniform> state: TrailState;

fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let x =  (p.x / screen.size.x) * 2.0 - 1.0;
    let y = -(p.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

struct VSOut {
    @builtin(position) pos: vec4<f32>,
    @location(0) color: vec4<f32>,
};

@vertex
fn vs_trail(
    @builtin(vertex_index) vi: u32,
    @builtin(instance_index) body: u32,
) -> VSOut {
    let seg = vi / 6u;
    let tri = vi % 6u;

    let i0 = (state.head + seg) % state.cap;
    let i1 = (state.head + seg + 1u) % state.cap;

    let idx0 = i0 * state.n_bodies + body;
    let idx1 = i1 * state.n_bodies + body;

    let p0 = positions[idx0];
    let p1 = positions[idx1];

    // descarta segmentos inválidos (NaN) — isNan não existe em WGSL; NaN != NaN é o teste correto
    if any(p0 != p0) || any(p1 != p1) {
        var out: VSOut;
        out.pos = vec4<f32>(0.0);
        out.color = vec4<f32>(0.0);
        return out;
    }

    let screen_p0 = state.center + p0 * state.scale;
    let screen_p1 = state.center + p1 * state.scale;

    let dir = screen_p1 - screen_p0;
    let len = max(length(dir), 1e-5);
    let tangent = dir / len;
    let normal = vec2<f32>(-tangent.y, tangent.x);

    let half_width = state.trail_width;

    var uv = array<vec2<f32>, 6>(
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0,  1.0),
    );

    let corner = uv[tri];
    let pos = screen_p0 + tangent * (len * corner.x) + normal * (half_width * corner.y);

    let base_color = colors[body];

    // fade ao longo do trail
    let t = f32(seg) / f32(state.cap);
    let alpha = base_color.a * (1.0 - t);

    var out: VSOut;
    out.pos = to_ndc(pos);
    out.color = vec4<f32>(base_color.rgb, alpha);
    return out;
}

@fragment
fn fs_trail(in: VSOut) -> @location(0) vec4<f32> {
    if in.color.a <= 0.001 {
        discard;
    }
    return in.color;
}
"#;
