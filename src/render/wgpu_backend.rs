use std::mem::size_of;
use std::sync::Arc;

use bytemuck::{Pod, Zeroable};
use egui;

use crate::render::TrailRenderer;
use crate::render::grid_renderer::GridRenderer;

const MIN_BUFFER_CAPACITY: u32 = 256;

// ── GPU data structures ───────────────────────────────────────────────────────

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ScreenUniform {
    /// Canvas dimensions in logical pixels (width, height).
    size: [f32; 2],
    /// Canvas origin in logical pixels (rect.min.x, rect.min.y).
    viewport_min: [f32; 2],
}

/// Global lighting state uploaded once per frame.
///
/// Used as a fallback when no per-instance light position is provided,
/// and to drive the ambient fill term for all bodies.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct LightingUniform {
    /// Screen-space position of the primary (brightest) light source in
    /// logical pixels. Used for bodies whose `light_pos` instance attribute
    /// is `[NaN, NaN]`.
    primary_light_pos: [f32; 2],
    /// Diffuse intensity of the primary light source [0, 1].
    intensity: f32,
    /// Ambient fill term [0, 1] — prevents fully unlit dark sides.
    ambient: f32,
}

impl Default for LightingUniform {
    fn default() -> Self {
        Self { primary_light_pos: [f32::NAN, f32::NAN], intensity: 0.6, ambient: 0.35 }
    }
}

/// Per-body circle instance.
///
/// `light_pos` carries the screen-space position of the nearest luminous
/// body to this instance. When set to `[NaN, NaN]` the fragment shader
/// falls back to the global `LightingUniform::primary_light_pos`.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CircleInstance {
    center: [f32; 2],
    outer_radius: f32,
    inner_radius: f32,
    /// Screen-space position of the nearest luminous body [px].
    /// `[NaN, NaN]` → use global primary light.
    light_pos: [f32; 2],
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

// ── GpuResources ─────────────────────────────────────────────────────────────

struct GpuResources {
    circle_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,

    screen_buf: wgpu::Buffer,
    lighting_buf: wgpu::Buffer,
    bind_group_screen: wgpu::BindGroup,
    bind_group_layout_screen: wgpu::BindGroupLayout,

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
            label: Some("screen_uniform"),
            size: size_of::<ScreenUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let lighting_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("lighting_uniform"),
            size: size_of::<LightingUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group_screen = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("screen_bg"),
            layout: &screen_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: screen_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: lighting_buf.as_entire_binding() },
            ],
        });

        let (circle_buf, circle_cap) =
            alloc_instance_buf::<CircleInstance>(device, MIN_BUFFER_CAPACITY, "circle");
        let (line_buf, line_cap) =
            alloc_instance_buf::<LineInstance>(device, MIN_BUFFER_CAPACITY, "line");

        Self {
            circle_pipeline,
            line_pipeline,
            screen_buf,
            lighting_buf,
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
        lighting: LightingUniform,
        circles: &[CircleInstance],
        lines: &[LineInstance],
    ) -> (u32, u32) {
        queue.write_buffer(&self.screen_buf, 0, bytemuck::bytes_of(&screen));
        queue.write_buffer(&self.lighting_buf, 0, bytemuck::bytes_of(&lighting));

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

// ── WgpuBackend ───────────────────────────────────────────────────────────────

/// Registered light source: screen-space position and luminosity weight.
///
/// The renderer uses this to assign per-body `light_pos` instance attributes
/// and to choose the global primary light for the `LightingUniform`.
#[derive(Clone, Copy)]
pub struct LightSource {
    /// Screen-space position in logical pixels.
    pub screen_pos: [f32; 2],
    /// Relative luminosity — used to find the brightest source and to
    /// weight multi-source blending in future extensions.
    pub luminosity: f32,
}

pub struct WgpuBackend {
    circles: Vec<CircleInstance>,
    lines: Vec<LineInstance>,

    /// Light sources registered this frame, derived from luminous bodies.
    /// Cleared at [`begin`] and repopulated by the app layer before draw.
    light_sources: Vec<LightSource>,

    /// Global lighting parameters (ambient, intensity).
    /// `primary_light_pos` is overwritten each frame from `light_sources`.
    lighting: LightingUniform,

    gpu: Option<GpuResources>,
    trail: Option<TrailRenderer>,
    grid: Option<GridRenderer>,

    pub trail_buffer: Option<Arc<crate::render::trail_buffer::TrailBuffer>>,
    pub trail_visibility: Option<Vec<bool>>,
    pub center: [f32; 2],
    pub scale: f32,
    pub show_grid: bool,
    pub trail_width: f32,
    pub trail_decay_k: f32,
    pub trail_tail_desaturate: f32,
}

impl WgpuBackend {
    pub fn new() -> Self {
        Self {
            circles: Vec::new(),
            lines: Vec::new(),
            light_sources: Vec::new(),
            lighting: LightingUniform::default(),
            gpu: None,
            trail: None,
            grid: None,

            trail_buffer: None,
            trail_visibility: None,
            center: [0.0, 0.0],
            scale: 1.0,
            show_grid: true,
            trail_width: 1.5,
            trail_decay_k: 6.0,
            trail_tail_desaturate: 0.5,
        }
    }

    pub fn begin(&mut self) {
        self.circles.clear();
        self.lines.clear();
        self.light_sources.clear();
    }

    /// Draw grid coordinate labels using egui on top of the GPU grid.
    ///
    /// Call this after the GPU paint callback, passing the same `rect` used for
    /// the wgpu viewport. `unit` is appended to x-axis tick labels (e.g. "AU").
    pub fn draw_labels_overlay(&self, ui: &egui::Ui, rect: egui::Rect, unit: &str) {
        if !self.show_grid {
            return;
        }
        if let Some(grid) = &self.grid {
            grid.draw_labels(ui, self.center, self.scale, rect, unit);
        }
    }

    // ── Lighting API ──────────────────────────────────────────────────────────

    /// Registers a luminous body as a light source for this frame.
    ///
    /// Call once per luminous body after [`begin`] and before [`draw_circle`].
    /// The backend uses these to:
    /// - select the global primary light (brightest source → `LightingUniform`)
    /// - assign per-instance `light_pos` to each circle (nearest source)
    pub fn add_light_source(&mut self, source: LightSource) {
        self.light_sources.push(source);
    }

    /// Sets the ambient fill and diffuse intensity for the global lighting model.
    ///
    /// Defaults: `ambient = 0.35`, `intensity = 0.6`.
    pub fn set_lighting_params(&mut self, ambient: f32, intensity: f32) {
        self.lighting.ambient = ambient.clamp(0.0, 1.0);
        self.lighting.intensity = intensity.clamp(0.0, 1.0);
    }

    /// Finds the nearest registered light source to `pos` in screen space.
    ///
    /// Returns `[NaN, NaN]` when no sources are registered, which causes the
    /// shader to fall back to the global `primary_light_pos`.
    fn nearest_light(&self, pos: [f32; 2]) -> [f32; 2] {
        self.light_sources
            .iter()
            .min_by(|a, b| {
                let da = dist2(pos, a.screen_pos);
                let db = dist2(pos, b.screen_pos);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|s| s.screen_pos)
            .unwrap_or([f32::NAN, f32::NAN])
    }

    /// Selects the primary light source (highest luminosity) for the global
    /// `LightingUniform`. Falls back to `[NaN, NaN]` when none are registered.
    fn primary_light_pos(&self) -> [f32; 2] {
        self.light_sources
            .iter()
            .max_by(|a, b| {
                a.luminosity.partial_cmp(&b.luminosity).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|s| s.screen_pos)
            .unwrap_or([f32::NAN, f32::NAN])
    }

    // ── Draw API ──────────────────────────────────────────────────────────────

    /// Submits a filled circle with per-instance lighting from the nearest
    /// registered [`LightSource`].
    pub fn draw_circle(&mut self, pos: [f32; 2], radius: f32, color: [u8; 3]) {
        let light_pos = self.nearest_light(pos);
        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: radius.max(0.5),
            inner_radius: 0.0,
            light_pos,
            color: rgba_u8_to_f32(color[0], color[1], color[2], 255),
        });
    }

    /// Submits an annular ring (stroke) with per-instance lighting.
    pub fn draw_circle_stroke(&mut self, pos: [f32; 2], radius: f32, width: f32, color: [u8; 4]) {
        let half = (width * 0.5).max(0.25);
        let outer = (radius + half).max(0.5);
        let inner = (radius - half).clamp(0.0, outer);
        let light_pos = self.nearest_light(pos);

        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: outer,
            inner_radius: inner,
            light_pos,
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

    // ── GPU setup ─────────────────────────────────────────────────────────────

    pub fn ensure_gpu(&mut self, device: &wgpu::Device, format: wgpu::TextureFormat) {
        if self.gpu.is_none() {
            self.gpu = Some(GpuResources::init(device, format));
        }
        if self.trail.is_none() {
            let gpu = self.gpu.as_ref().unwrap();
            self.trail = Some(TrailRenderer::new(device, &gpu.bind_group_layout_screen, format));
        }
        if self.grid.is_none() {
            self.grid = Some(GridRenderer::new(device, format));
        }
    }

    // ── Render ────────────────────────────────────────────────────────────────

    pub fn render_frame(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        pass: &mut wgpu::RenderPass<'_>,
        screen: [f32; 2],
        viewport_min: [f32; 2],
        format: wgpu::TextureFormat,
        center: [f32; 2],
        scale: f32,
    ) {
        self.ensure_gpu(device, format);
        self.center = center;
        self.scale = scale;

        // Resolve primary light for the global uniform.
        let mut lighting = self.lighting;
        lighting.primary_light_pos = self.primary_light_pos();

        let screen_uniform = ScreenUniform { size: screen, viewport_min };

        let (circle_count, line_count) = {
            let gpu = self.gpu.as_mut().unwrap();
            gpu.upload(device, queue, screen_uniform, lighting, &self.circles, &self.lines)
        };

        if self.show_grid {
            if let Some(grid) = &self.grid {
                grid.upload(queue, center, scale, screen);
                grid.draw(pass);
            }
        }

        if let (Some(trail_renderer), Some(trail_buf)) =
            (self.trail.as_mut(), self.trail_buffer.as_deref())
        {
            let gpu = self.gpu.as_ref().unwrap();
            trail_renderer.upload(
                device,
                queue,
                trail_buf,
                self.trail_visibility.as_deref(),
                center,
                scale,
                self.trail_width,
                self.trail_decay_k,
                self.trail_tail_desaturate,
            );
            trail_renderer.draw(pass, &gpu.bind_group_screen);
        }

        {
            let gpu = self.gpu.as_ref().unwrap();
            gpu.draw(pass, circle_count, line_count);
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

#[inline]
fn dist2(a: [f32; 2], b: [f32; 2]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy
}

#[inline]
fn rgba_u8_to_f32(r: u8, g: u8, b: u8, a: u8) -> [f32; 4] {
    [r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0, a as f32 / 255.0]
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
    let (new_buf, new_c) = alloc_instance_buf::<T>(device, new_cap, label);

    if *cap > 0 {
        let copy_bytes = *cap as u64 * size_of::<T>() as u64;
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("instance_buf_grow"),
        });
        enc.copy_buffer_to_buffer(buf, 0, &new_buf, 0, copy_bytes);
        queue.submit([enc.finish()]);
    }

    *buf = new_buf;
    *cap = new_c;
}

// ── Pipeline builders ─────────────────────────────────────────────────────────

fn build_bind_group_layout(device: &wgpu::Device) -> wgpu::BindGroupLayout {
    device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("grav_sim::bgl"),
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::FRAGMENT,
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

fn build_shader(device: &wgpu::Device) -> wgpu::ShaderModule {
    device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("grav_sim::primitives"),
        source: wgpu::ShaderSource::Wgsl(PRIMITIVES_SHADER.into()),
    })
}

fn build_circle_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    layout: &wgpu::PipelineLayout,
    format: wgpu::TextureFormat,
) -> wgpu::RenderPipeline {
    // center(2) + outer_r(1) + inner_r(1) + light_pos(2) + _pad removed + color(4)
    let attrs = wgpu::vertex_attr_array![
        0 => Float32x2,   // center
        1 => Float32,     // outer_radius
        2 => Float32,     // inner_radius
        3 => Float32x2,   // light_pos
        4 => Float32x4    // color
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

const PRIMITIVES_SHADER: &str = r#"

struct ScreenUniform {
    size:         vec2<f32>,
    viewport_min: vec2<f32>,
};

/// Global lighting fallback and ambient parameters.
struct LightingUniform {
    /// Screen-space position of the primary light source.
    /// Components are NaN when no source is registered.
    primary_light_pos: vec2<f32>,
    /// Diffuse intensity [0, 1].
    intensity:         f32,
    /// Ambient fill [0, 1].
    ambient:           f32,
};

@group(0) @binding(0) var<uniform> screen:   ScreenUniform;
@group(0) @binding(1) var<uniform> lighting: LightingUniform;

fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let local = p - screen.viewport_min;
    let x =  (local.x / screen.size.x) * 2.0 - 1.0;
    let y = -(local.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

// ── CIRCLES ──────────────────────────────────────────────────────────────────

struct CircleVarying {
    @builtin(position) clip_pos:    vec4<f32>,
    /// Fragment position in logical pixels — used to compute light direction.
    @location(0)       screen_pos:  vec2<f32>,
    @location(1)       local:       vec2<f32>,
    @location(2)       inner_ratio: f32,
    /// Per-instance nearest light position in screen space.
    /// NaN signals "use global primary".
    @location(3)       light_pos:   vec2<f32>,
    @location(4)       color:       vec4<f32>,
};

@vertex
fn vs_circle(
    @builtin(vertex_index) vi:           u32,
    @location(0)           center:       vec2<f32>,
    @location(1)           outer_radius: f32,
    @location(2)           inner_radius: f32,
    @location(3)           light_pos:    vec2<f32>,
    @location(4)           color:        vec4<f32>,
) -> CircleVarying {
    var quad = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );

    let local      = quad[vi];
    let world_pos  = center + local * outer_radius;

    var out: CircleVarying;
    out.clip_pos    = to_ndc(world_pos);
    out.screen_pos  = world_pos;
    out.local       = local;
    out.inner_ratio = select(0.0, inner_radius / outer_radius, outer_radius > 0.0);
    out.light_pos   = light_pos;
    out.color       = color;
    return out;
}

@fragment
fn fs_circle(in: CircleVarying) -> @location(0) vec4<f32> {
    let r  = length(in.local);
    let aa = fwidth(r);

    let outer      = 1.0 - smoothstep(1.0 - aa, 1.0 + aa, r);
    let inner_r    = in.inner_ratio;
    let inner_mask = select(1.0, smoothstep(inner_r - aa, inner_r + aa, r), inner_r > 0.0);

    var alpha = in.color.a * outer * inner_mask;
    if alpha <= 0.001 { discard; }

    // ── Light source resolution ───────────────────────────────────────────
    // Use the per-instance light_pos when valid; fall back to the global
    // primary when either component is NaN (no per-instance source).
    var lp = in.light_pos;
    if (lp.x != lp.x) {   // NaN check: NaN != NaN is true in WGSL
        lp = lighting.primary_light_pos;
    }

    // ── Lighting model ────────────────────────────────────────────────────
    // When no light source at all is available (lp still NaN after fallback),
    // skip directional lighting and use ambient only.
    var diffuse = 0.0;
    if (lp.x == lp.x) {
        // Direction from fragment surface point to the light source.
        let to_light = normalize(lp - in.screen_pos);
        // 2D surface normal: the local coord of the fragment on the unit disc.
        let n        = normalize(in.local);
        diffuse      = clamp(dot(n, to_light), 0.0, 1.0);
    }

    let light_f   = lighting.ambient
                  + (1.0 - lighting.ambient) * lighting.intensity * diffuse;

    // Subtle limb darkening — depth cue independent of light direction.
    let edge      = smoothstep(0.7, 1.0, r);
    let edge_dark = mix(1.0, 0.75, edge);

    let lit_color = in.color.rgb * light_f * edge_dark;
    return vec4<f32>(lit_color, alpha);
}

// ── LINES ─────────────────────────────────────────────────────────────────────

struct LineVarying {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0)       uv:       vec2<f32>,
    @location(1)       color:    vec4<f32>,
};

@vertex
fn vs_line(
    @builtin(vertex_index) vi:         u32,
    @location(0)           p0:         vec2<f32>,
    @location(1)           p1:         vec2<f32>,
    @location(2)           half_width: f32,
    @location(3)           _pad:       f32,
    @location(4)           color:      vec4<f32>,
) -> LineVarying {
    let dir     = p1 - p0;
    let len     = max(length(dir), 1e-5);
    let tangent = dir / len;
    let normal  = vec2<f32>(-tangent.y, tangent.x);

    var uv = array<vec2<f32>, 6>(
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0,  1.0),
    );

    let corner = uv[vi];
    let pos    = p0 + tangent * (len * corner.x) + normal * (half_width * corner.y);

    var out: LineVarying;
    out.clip_pos = to_ndc(pos);
    out.uv       = corner;
    out.color    = color;
    return out;
}

@fragment
fn fs_line(in: LineVarying) -> @location(0) vec4<f32> {
    let edge  = abs(in.uv.y);
    let aa    = fwidth(edge);
    let alpha = 1.0 - smoothstep(1.0 - aa, 1.0 + aa, edge);

    let final_alpha = in.color.a * alpha;
    if final_alpha <= 0.001 { discard; }

    return vec4<f32>(in.color.rgb, final_alpha);
}
"#;
