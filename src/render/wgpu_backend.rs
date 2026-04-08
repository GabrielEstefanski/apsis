use std::sync::Mutex;

use bytemuck::{Pod, Zeroable};
use eframe::egui::{self, Rect};
use eframe::egui_wgpu::{self, CallbackTrait};
use wgpu::util::DeviceExt;

use crate::render::RenderBackend;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ScreenUniform {
    size: [f32; 2],
    _pad: [f32; 2],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CircleInstance {
    center: [f32; 2],
    outer_radius: f32,
    inner_radius: f32,
    _pad: f32,
    color: [f32; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct LineInstance {
    from: [f32; 2],
    to: [f32; 2],
    half_width: f32,
    _pad: f32,
    color: [f32; 4],
}

struct SharedGpuResources {
    bind_group_layout: wgpu::BindGroupLayout,
    circle_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,
}

struct PreparedFrame {
    bind_group: wgpu::BindGroup,
    circle_buffer: Option<wgpu::Buffer>,
    circle_count: u32,
    line_buffer: Option<wgpu::Buffer>,
    line_count: u32,
}

struct WgpuPrimitivesCallback {
    circles: Vec<CircleInstance>,
    lines: Vec<LineInstance>,
    prepared: Mutex<Option<PreparedFrame>>,
}

fn make_shared_resources(device: &wgpu::Device) -> SharedGpuResources {
    let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("gravity_sim_wgpu_bind_group_layout"),
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
    });

    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("gravity_sim_wgpu_primitives_shader"),
        source: wgpu::ShaderSource::Wgsl(PRIMITIVES_SHADER.into()),
    });

    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("gravity_sim_wgpu_pipeline_layout"),
        bind_group_layouts: &[Some(&bind_group_layout)],
        immediate_size: 0,
    });

    let circle_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("gravity_sim_wgpu_circle_pipeline"),
        layout: Some(&pipeline_layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: Some("vs_circle"),
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<CircleInstance>() as u64,
                step_mode: wgpu::VertexStepMode::Instance,
                attributes: &wgpu::vertex_attr_array![
                    0 => Float32x2,
                    1 => Float32,
                    2 => Float32,
                    3 => Float32,
                    4 => Float32x4
                ],
            }],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: Some("fs_circle"),
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
    });

    let line_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("gravity_sim_wgpu_line_pipeline"),
        layout: Some(&pipeline_layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: Some("vs_line"),
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<LineInstance>() as u64,
                step_mode: wgpu::VertexStepMode::Instance,
                attributes: &wgpu::vertex_attr_array![
                    0 => Float32x2,
                    1 => Float32x2,
                    2 => Float32,
                    3 => Float32,
                    4 => Float32x4
                ],
            }],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: Some("fs_line"),
            targets: &[Some(wgpu::ColorTargetState {
                format: wgpu::TextureFormat::Bgra8Unorm,
                blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                write_mask: wgpu::ColorWrites::ALL,
            })],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::LineList,
            ..Default::default()
        },
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        multiview_mask: None,
        cache: None,
    });

    SharedGpuResources {
        bind_group_layout,
        circle_pipeline,
        line_pipeline,
    }
}

impl CallbackTrait for WgpuPrimitivesCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        _queue: &wgpu::Queue,
        screen_descriptor: &egui_wgpu::ScreenDescriptor,
        _egui_encoder: &mut wgpu::CommandEncoder,
        callback_resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        if callback_resources.get::<SharedGpuResources>().is_none() {
            callback_resources.insert(make_shared_resources(device));
        }
        let shared = callback_resources
            .get::<SharedGpuResources>()
            .expect("wgpu shared resources missing");

        let screen_uniform = ScreenUniform {
            size: [
                screen_descriptor.size_in_pixels[0] as f32 / screen_descriptor.pixels_per_point,
                screen_descriptor.size_in_pixels[1] as f32 / screen_descriptor.pixels_per_point,
            ],
            _pad: [0.0, 0.0],
        };
        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("gravity_sim_wgpu_uniform"),
            contents: bytemuck::bytes_of(&screen_uniform),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("gravity_sim_wgpu_bind_group"),
            layout: &shared.bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let circle_buffer = if self.circles.is_empty() {
            None
        } else {
            Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("gravity_sim_wgpu_circle_buffer"),
                contents: bytemuck::cast_slice(&self.circles),
                usage: wgpu::BufferUsages::VERTEX,
            }))
        };

        let line_buffer = if self.lines.is_empty() {
            None
        } else {
            Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("gravity_sim_wgpu_line_buffer"),
                contents: bytemuck::cast_slice(&self.lines),
                usage: wgpu::BufferUsages::VERTEX,
            }))
        };

        *self.prepared.lock().expect("wgpu prepared lock poisoned") = Some(PreparedFrame {
            bind_group,
            circle_buffer,
            circle_count: self.circles.len() as u32,
            line_buffer,
            line_count: self.lines.len() as u32,
        });

        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        callback_resources: &egui_wgpu::CallbackResources,
    ) {
        let guard = self.prepared.lock().expect("wgpu prepared lock poisoned");
        let Some(prepared) = guard.as_ref() else {
            return;
        };
        let Some(shared) = callback_resources.get::<SharedGpuResources>() else {
            return;
        };

        if let Some(line_buffer) = prepared.line_buffer.as_ref() {
            render_pass.set_pipeline(&shared.line_pipeline);
            render_pass.set_bind_group(0, &prepared.bind_group, &[]);
            render_pass.set_vertex_buffer(0, line_buffer.slice(..));
            render_pass.draw(0..6, 0..prepared.line_count);
        }

        if let Some(circle_buffer) = prepared.circle_buffer.as_ref() {
            render_pass.set_pipeline(&shared.circle_pipeline);
            render_pass.set_bind_group(0, &prepared.bind_group, &[]);
            render_pass.set_vertex_buffer(0, circle_buffer.slice(..));
            render_pass.draw(0..6, 0..prepared.circle_count);
        }
    }
}

pub struct WgpuBackend<'a> {
    ui: &'a egui::Ui,
    rect: Rect,
    circles: Vec<CircleInstance>,
    lines: Vec<LineInstance>,
}

impl<'a> WgpuBackend<'a> {
    pub fn new(ui: &'a egui::Ui, rect: Rect) -> Self {
        Self {
            ui,
            rect,
            circles: Vec::new(),
            lines: Vec::new(),
        }
    }
}

impl RenderBackend for WgpuBackend<'_> {
    fn begin(&mut self) {}

    fn draw_circle(&mut self, pos: [f32; 2], radius: f32, color: [u8; 3]) {
        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: radius.max(0.5),
            inner_radius: 0.0,
            _pad: 0.0,
            color: [
                color[0] as f32 / 255.0,
                color[1] as f32 / 255.0,
                color[2] as f32 / 255.0,
                1.0,
            ],
        });
    }

    fn draw_circle_stroke(&mut self, pos: [f32; 2], radius: f32, width: f32, color: [u8; 4]) {
        let half = (width * 0.5).max(0.25);
        let outer = (radius + half).max(0.5);
        let inner = (radius - half).max(0.0);
        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: outer,
            inner_radius: inner.min(outer),
            _pad: 0.0,
            color: [
                color[0] as f32 / 255.0,
                color[1] as f32 / 255.0,
                color[2] as f32 / 255.0,
                color[3] as f32 / 255.0,
            ],
        });
    }

    fn draw_line_segment(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: [u8; 4]) {
        let color = [
            color[0] as f32 / 255.0,
            color[1] as f32 / 255.0,
            color[2] as f32 / 255.0,
            color[3] as f32 / 255.0,
        ];
        self.lines.push(LineInstance {
            from,
            to,
            half_width: (width * 0.5).max(0.25),
            _pad: 0.0,
            color,
        });
    }

    fn end(&mut self) {
        let callback = egui_wgpu::Callback::new_paint_callback(
            self.rect,
            WgpuPrimitivesCallback {
                circles: std::mem::take(&mut self.circles),
                lines: std::mem::take(&mut self.lines),
                prepared: Mutex::new(None),
            },
        );

        self.ui.painter().add(callback);
    }
}

const PRIMITIVES_SHADER: &str = r#"
struct ScreenUniform {
    size: vec2<f32>,
    _pad: vec2<f32>,
};

@group(0) @binding(0)
var<uniform> screen: ScreenUniform;

struct CircleOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) local: vec2<f32>,
    @location(1) inner_ratio: f32,
    @location(2) color: vec4<f32>,
};

fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let x = (p.x / screen.size.x) * 2.0 - 1.0;
    let y = 1.0 - (p.y / screen.size.y) * 2.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

@vertex
fn vs_circle(
    @builtin(vertex_index) vertex_index: u32,
    @location(0) center: vec2<f32>,
    @location(1) outer_radius: f32,
    @location(2) inner_radius: f32,
    @location(3) _pad: f32,
    @location(4) color: vec4<f32>,
) -> CircleOut {
    var quad = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );
    let local = quad[vertex_index];
    let pos = center + local * outer_radius;

    var out: CircleOut;
    out.clip_pos = to_ndc(pos);
    out.local = local;
    out.inner_ratio = select(0.0, inner_radius / outer_radius, outer_radius > 0.0);
    out.color = color;
    return out;
}

@fragment
fn fs_circle(in: CircleOut) -> @location(0) vec4<f32> {
    let r2 = dot(in.local, in.local);
    if (r2 > 1.0 || r2 < in.inner_ratio * in.inner_ratio) {
        discard;
    }
    return in.color;
}

struct LineOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) color: vec4<f32>,
};

@vertex
fn vs_line(
    @builtin(vertex_index) vertex_index: u32,
    @location(0) p0: vec2<f32>,
    @location(1) p1: vec2<f32>,
    @location(2) half_width: f32,
    @location(3) _pad: f32,
    @location(4) color: vec4<f32>,
) -> LineOut {
    let dir = p1 - p0;
    let len = max(length(dir), 1e-5);
    let tangent = dir / len;
    let normal = vec2<f32>(-tangent.y, tangent.x) * half_width;

    var quad = array<vec2<f32>, 6>(
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0,  1.0),
    );

    let uv = quad[vertex_index];
    let along = p0 + dir * uv.x;
    let pos = along + normal * uv.y;

    var out: LineOut;
    out.clip_pos = to_ndc(pos);
    out.color = color;
    return out;
}

@fragment
fn fs_line(in: LineOut) -> @location(0) vec4<f32> {
    return in.color;
}
"#;
