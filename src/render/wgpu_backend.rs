use std::mem::size_of;

use bytemuck::{Pod, Zeroable};

use crate::render::TrailRenderer;
use crate::render::grid_renderer::GridRenderer;

const MIN_BUFFER_CAPACITY: u32 = 256;

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

struct GpuResources {
    // pipelines
    circle_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,

    // screen uniform
    screen_buf: wgpu::Buffer,
    bind_group_screen: wgpu::BindGroup,
    bind_group_layout_screen: wgpu::BindGroupLayout,

    // instance buffers
    circle_buf: wgpu::Buffer,
    circle_cap: u32,

    line_buf: wgpu::Buffer,
    line_cap: u32,
}

impl GpuResources {
    fn init(device: &wgpu::Device, format: wgpu::TextureFormat) -> Self {
        let screen_bgl = build_bind_group_layout(device);

        let shader = build_shader(device);

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pipeline_layout"),
            bind_group_layouts: &[Some(&screen_bgl)],
            immediate_size: 0,
        });

        let circle_pipeline = build_circle_pipeline(device, &shader, &pipeline_layout, format);

        let line_pipeline = build_line_pipeline(device, &shader, &pipeline_layout, format);

        let screen_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("screen"),
            size: size_of::<ScreenUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group_screen = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("screen_bg"),
            layout: &screen_bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: screen_buf.as_entire_binding(),
            }],
        });

        let (circle_buf, circle_cap) =
            alloc_instance_buf::<CircleInstance>(device, MIN_BUFFER_CAPACITY, "circle");

        let (line_buf, line_cap) =
            alloc_instance_buf::<LineInstance>(device, MIN_BUFFER_CAPACITY, "line");

        Self {
            circle_pipeline,
            line_pipeline,
            screen_buf,
            bind_group_screen,
            bind_group_layout_screen: screen_bgl,
            circle_buf,
            circle_cap,
            line_buf,
            line_cap,
        }
    }

    fn upload(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        screen: ScreenUniform,
        circles: &[CircleInstance],
        lines: &[LineInstance],
    ) -> (u32, u32) {
        queue.write_buffer(&self.screen_buf, 0, bytemuck::bytes_of(&screen));

        if !circles.is_empty() {
            ensure_instance_capacity::<CircleInstance>(
                device,
                queue,
                &mut self.circle_buf,
                &mut self.circle_cap,
                circles.len() as u32,
                "circle",
            );
            queue.write_buffer(&self.circle_buf, 0, bytemuck::cast_slice(circles));
        }

        if !lines.is_empty() {
            ensure_instance_capacity::<LineInstance>(
                device,
                queue,
                &mut self.line_buf,
                &mut self.line_cap,
                lines.len() as u32,
                "line",
            );
            queue.write_buffer(&self.line_buf, 0, bytemuck::cast_slice(lines));
        }

        (circles.len() as u32, lines.len() as u32)
    }

    fn draw(&self, pass: &mut wgpu::RenderPass<'_>, circle_count: u32, line_count: u32) {
        if line_count > 0 {
            pass.set_pipeline(&self.line_pipeline);
            pass.set_bind_group(0, &self.bind_group_screen, &[]);
            pass.set_vertex_buffer(0, self.line_buf.slice(..));
            pass.draw(0..6, 0..line_count);
        }

        if circle_count > 0 {
            pass.set_pipeline(&self.circle_pipeline);
            pass.set_bind_group(0, &self.bind_group_screen, &[]);
            pass.set_vertex_buffer(0, self.circle_buf.slice(..));
            pass.draw(0..6, 0..circle_count);
        }
    }
}

pub struct WgpuBackend {
    circles: Vec<CircleInstance>,
    lines: Vec<LineInstance>,

    gpu: Option<GpuResources>,
    trail: Option<TrailRenderer>,
    grid: Option<GridRenderer>,

    pub trail_buffer: Option<crate::core::trail_buffer::TrailBuffer>,
    pub center: [f32; 2],
    pub scale: f32,
    pub show_grid: bool,
    pub trail_width: f32,
}

impl WgpuBackend {
    pub fn new() -> Self {
        Self {
            circles: Vec::new(),
            lines: Vec::new(),
            gpu: None,
            trail: None,
            grid: None,

            trail_buffer: None,
            center: [0.0, 0.0],
            scale: 1.0,
            show_grid: true,
            trail_width: 1.5,
        }
    }

    pub fn begin(&mut self) {
        self.circles.clear();
        self.lines.clear();
    }

    // ── DRAW API ────────────────────────────────────────────

    pub fn draw_circle(&mut self, pos: [f32; 2], radius: f32, color: [u8; 3]) {
        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: radius.max(0.5),
            inner_radius: 0.0,
            _pad: 0.0,
            color: rgba_u8_to_f32(color[0], color[1], color[2], 255),
        });
    }

    pub fn draw_circle_stroke(&mut self, pos: [f32; 2], radius: f32, width: f32, color: [u8; 4]) {
        let half = (width * 0.5).max(0.25);
        let outer = (radius + half).max(0.5);
        let inner = (radius - half).clamp(0.0, outer);

        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: outer,
            inner_radius: inner,
            _pad: 0.0,
            color: rgba_u8_to_f32(color[0], color[1], color[2], color[3]),
        });
    }

    pub fn draw_line_segment(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: [u8; 4]) {
        self.lines.push(LineInstance {
            from,
            to,
            half_width: (width * 0.5).max(0.25),
            _pad: 0.0,
            color: rgba_u8_to_f32(color[0], color[1], color[2], color[3]),
        });
    }

    // ── GPU SETUP ───────────────────────────────────────────

    pub fn ensure_gpu(&mut self, device: &wgpu::Device, format: wgpu::TextureFormat) {
        if self.gpu.is_none() {
            self.gpu = Some(GpuResources::init(device, format));
        }

        if self.trail.is_none() {
            let gpu = self.gpu.as_ref().unwrap();
            self.trail = Some(TrailRenderer::new(
                device,
                &gpu.bind_group_layout_screen,
                format,
            ));
        }

        if self.grid.is_none() {
            self.grid = Some(GridRenderer::new(device, format));
        }
    }

    // ── RENDER ──────────────────────────────────────────────

    pub fn render_frame(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        pass: &mut wgpu::RenderPass<'_>,
        screen: [f32; 2],
        format: wgpu::TextureFormat,
        trail_buf: Option<&crate::core::trail_buffer::TrailBuffer>,
        center: [f32; 2],
        scale: f32,
    ) {
        self.ensure_gpu(device, format);

        self.center = center;
        self.scale = scale;

        // ── 1. Grid (fundo) ─────────────────────────────────
        if self.show_grid {
            if let Some(grid) = &self.grid {
                grid.upload(queue, center, scale, screen);
                grid.draw(pass);
            }
        }

        let (circle_count, line_count) = {
            let gpu = self.gpu.as_mut().unwrap();
            let screen_uniform = ScreenUniform {
                size: screen,
                _pad: [0.0; 2],
            };
            gpu.upload(device, queue, screen_uniform, &self.circles, &self.lines)
        };

        let gpu = self.gpu.as_ref().unwrap();

        // ── 2. Trails ────────────────────────────────────────
        if let (Some(trail), Some(buf)) = (self.trail.as_mut(), trail_buf) {
            trail.upload(device, queue, buf, center, scale, self.trail_width);
            trail.draw(pass, &gpu.bind_group_screen);
        }

        // ── 3. Corpos (circles / lines) ──────────────────────
        gpu.draw(pass, circle_count, line_count);
    }
}

#[inline]
fn rgba_u8_to_f32(r: u8, g: u8, b: u8, a: u8) -> [f32; 4] {
    [
        r as f32 / 255.0,
        g as f32 / 255.0,
        b as f32 / 255.0,
        a as f32 / 255.0,
    ]
}

fn alloc_instance_buf<T: Pod>(
    device: &wgpu::Device,
    capacity: u32,
    label: &str,
) -> (wgpu::Buffer, u32) {
    let buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some(label),
        size: capacity as u64 * size_of::<T>() as u64,
        usage: wgpu::BufferUsages::VERTEX
            | wgpu::BufferUsages::COPY_DST
            | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    (buf, capacity)
}

fn ensure_instance_capacity<T: Pod>(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    buf: &mut wgpu::Buffer,
    cap: &mut u32,
    needed: u32,
    label: &str,
) {
    if needed <= *cap {
        return;
    }

    let new_cap = needed.next_power_of_two().max(MIN_BUFFER_CAPACITY);
    let (new_buf, new_cap) = alloc_instance_buf::<T>(device, new_cap, label);

    // Copy existing GPU data into the new larger buffer so growth is seamless.
    // Instance buffers are fully rewritten each frame via write_buffer, so this
    // is mostly for correctness; the real win is avoiding the old buffer's
    // data being stale during the brief window between resize and upload.
    if *cap > 0 {
        let copy_bytes = *cap as u64 * size_of::<T>() as u64;
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("instance_buf_grow"),
        });
        enc.copy_buffer_to_buffer(buf, 0, &new_buf, 0, copy_bytes);
        queue.submit([enc.finish()]);
    }

    *buf = new_buf;
    *cap = new_cap;
}

// ── Pipeline / shader builders ────────────────────────────────────────────────

fn build_bind_group_layout(device: &wgpu::Device) -> wgpu::BindGroupLayout {
    device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("grav_sim::bgl"),
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

fn build_shader(device: &wgpu::Device) -> wgpu::ShaderModule {
    device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("grav_sim::primitives"),
        source: wgpu::ShaderSource::Wgsl(PRIMITIVES_SHADER.into()),
    })
}

fn build_pipeline_layout(
    device: &wgpu::Device,
    bgl: &wgpu::BindGroupLayout,
) -> wgpu::PipelineLayout {
    device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("grav_sim::pipeline_layout"),
        bind_group_layouts: &[Some(bgl)],
        immediate_size: 0,
    })
}

fn build_circle_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    layout: &wgpu::PipelineLayout,
    format: wgpu::TextureFormat,
) -> wgpu::RenderPipeline {
    let attrs = wgpu::vertex_attr_array![
        0 => Float32x2,
        1 => Float32,
        2 => Float32,
        3 => Float32,
        4 => Float32x4
    ];

    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("grav_sim::circle_pipeline"),
        layout: Some(layout),
        vertex: wgpu::VertexState {
            module: shader,
            entry_point: Some("vs_circle"),
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: size_of::<CircleInstance>() as u64,
                step_mode: wgpu::VertexStepMode::Instance,
                attributes: &attrs,
            }],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: shader,
            entry_point: Some("fs_circle"),
            targets: &[Some(wgpu::ColorTargetState {
                format,
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

fn build_line_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    layout: &wgpu::PipelineLayout,
    format: wgpu::TextureFormat,
) -> wgpu::RenderPipeline {
    let attrs = wgpu::vertex_attr_array![
        0 => Float32x2,
        1 => Float32x2,
        2 => Float32,
        3 => Float32,
        4 => Float32x4
    ];

    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("grav_sim::line_pipeline"),
        layout: Some(layout),
        vertex: wgpu::VertexState {
            module: shader,
            entry_point: Some("vs_line"),
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: size_of::<LineInstance>() as u64,
                step_mode: wgpu::VertexStepMode::Instance,
                attributes: &attrs,
            }],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: shader,
            entry_point: Some("fs_line"),
            targets: &[Some(wgpu::ColorTargetState {
                format,
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

/// Shared WGSL shader for circles and line segments.
///
/// ## Circle pass
///
/// Each instance generates a quad bounding the disc.  The fragment shader
/// discards pixels outside `[inner_ratio, 1.0]` in normalised disc space,
/// producing a solid disc or annular ring.
///
/// ## Line pass
///
/// Each instance generates a screen-aligned quad extruded along the segment
/// normal by `half_width`.  No anti-aliasing is applied; soft edges can be
/// added in a future pass.
const PRIMITIVES_SHADER: &str = r#"
struct ScreenUniform {
    size: vec2<f32>,
    _pad: vec2<f32>,
};

@group(0) @binding(0)
var<uniform> screen: ScreenUniform;

fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let x =  (p.x / screen.size.x) * 2.0 - 1.0;
    let y = -(p.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

// ── CIRCLES ─────────────────────────────────────────────────────────────── //

struct CircleVarying {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) local: vec2<f32>,
    @location(1) inner_ratio: f32,
    @location(2) color: vec4<f32>,
};

@vertex
fn vs_circle(
    @builtin(vertex_index) vi: u32,
    @location(0) center: vec2<f32>,
    @location(1) outer_radius: f32,
    @location(2) inner_radius: f32,
    @location(3) _pad: f32,
    @location(4) color: vec4<f32>,
) -> CircleVarying {
    var quad = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );

    let local = quad[vi];
    let pos   = center + local * outer_radius;

    var out: CircleVarying;
    out.clip_pos    = to_ndc(pos);
    out.local       = local;
    out.inner_ratio = select(0.0, inner_radius / outer_radius, outer_radius > 0.0);
    out.color       = color;
    return out;
}

@fragment
fn fs_circle(in: CircleVarying) -> @location(0) vec4<f32> {
    let r = length(in.local);
    let aa = fwidth(r);

    // Outer edge AA
    let outer = 1.0 - smoothstep(1.0 - aa, 1.0 + aa, r);

    // Inner cut (rings)
    let inner = in.inner_ratio;
    let inner_mask = select(1.0, smoothstep(inner - aa, inner + aa, r), inner > 0.0);

    var alpha = in.color.a * outer * inner_mask;
    if alpha <= 0.001 { discard; }

    // ── SHADING (melhora MUITO a qualidade visual) ──
    let light_dir = normalize(vec2<f32>(-0.6, -0.8));
    let n = normalize(in.local);
    let diffuse = clamp(dot(n, -light_dir), 0.0, 1.0);

    // Fake spherical lighting
    let lighting = 0.4 + 0.6 * diffuse;

    // Subtle edge darkening (depth cue)
    let edge = smoothstep(0.7, 1.0, r);
    let edge_dark = mix(1.0, 0.75, edge);

    let color = in.color.rgb * lighting * edge_dark;

    return vec4<f32>(color, alpha);
}

// ── LINES ─────────────────────────────────────────────────────────────── //

struct LineVarying {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) uv: vec2<f32>,
    @location(1) color: vec4<f32>,
};

@vertex
fn vs_line(
    @builtin(vertex_index) vi: u32,
    @location(0) p0: vec2<f32>,
    @location(1) p1: vec2<f32>,
    @location(2) half_width: f32,
    @location(3) _pad: f32,
    @location(4) color: vec4<f32>,
) -> LineVarying {
    let dir = p1 - p0;
    let len = max(length(dir), 1e-5);
    let tangent = dir / len;
    let normal = vec2<f32>(-tangent.y, tangent.x);

    var uv = array<vec2<f32>, 6>(
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0,  1.0),
    );

    let corner = uv[vi];
    let pos = p0 + tangent * (len * corner.x) + normal * (half_width * corner.y);

    var out: LineVarying;
    out.clip_pos = to_ndc(pos);
    out.uv = corner;
    out.color = color;
    return out;
}

@fragment
fn fs_line(in: LineVarying) -> @location(0) vec4<f32> {
    // Soft edge (AA)
    let edge = abs(in.uv.y);
    let aa = fwidth(edge);
    let alpha = 1.0 - smoothstep(1.0 - aa, 1.0 + aa, edge);

    let final_alpha = in.color.a * alpha;
    if final_alpha <= 0.001 { discard; }

    return vec4<f32>(in.color.rgb, final_alpha);
}
"#;
