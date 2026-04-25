use std::mem::size_of;
use std::sync::Arc;

use bytemuck::{Pod, Zeroable};
use egui;

use std::time::Instant;

use crate::render::TrailRenderer;
use crate::render::exposure::{ExposureState, decode_reduced_texel};
use crate::render::grid_renderer::GridRenderer;
use crate::render::hdr::{HDR_FORMAT, HdrTarget};
use crate::render::lighting::{LightingUniform, SceneLighting};
use crate::render::luminance_reducer::LuminanceReducer;
use crate::render::tonemap::TonemapPipeline;

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

/// Flat circle / ring primitive — used strictly for annotations (apsides
/// markers, selection rings, etc.) that must *not* be affected by scene
/// lighting. Bodies go through the dedicated [`BodyInstance`] pipeline.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CircleInstance {
    center: [f32; 2],
    outer_radius: f32,
    inner_radius: f32,
    color: [f32; 4],
}

/// Per-body sphere-from-disc instance.
///
/// Geometry data (`center_screen`, `radius_screen`) drives the 2D quad
/// expansion in the vertex stage. Lighting data (`center_world`, `albedo`,
/// `emissive`) lives alongside so the fragment shader can reconstruct a
/// sphere normal and compute a world-space Lambert term in one pass.
///
/// `albedo` is the diffuse base colour (RGBA, alpha controls fade for UI
/// blending). `emissive` is added unattenuated on top — stars set this to
/// their body colour and albedo near zero so their disc glows self-lit.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct BodyInstance {
    center_screen: [f32; 2],
    center_world: [f32; 3],
    radius_screen: f32,
    albedo: [f32; 4],
    emissive: [f32; 4],
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
    body_pipeline: wgpu::RenderPipeline,
    circle_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,

    screen_buf: wgpu::Buffer,
    lighting_buf: wgpu::Buffer,
    bind_group_screen: wgpu::BindGroup,
    bind_group_layout_screen: wgpu::BindGroupLayout,

    body_buf: wgpu::Buffer,
    body_cap: u32,
    circle_buf: wgpu::Buffer,
    circle_cap: u32,
    line_buf: wgpu::Buffer,
    line_cap: u32,
}

impl GpuResources {
    fn init(device: &wgpu::Device) -> Self {
        let screen_bgl = build_bind_group_layout(device);
        let shader = build_shader(device);

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pipeline_layout"),
            bind_group_layouts: &[Some(&screen_bgl)],
            immediate_size: 0,
        });

        // Scene primitives always render into the HDR offscreen target; the
        // swapchain format is handled by the tonemap pass.
        let body_pipeline = build_body_pipeline(device, &shader, &pipeline_layout, HDR_FORMAT);
        let circle_pipeline = build_circle_pipeline(device, &shader, &pipeline_layout, HDR_FORMAT);
        let line_pipeline = build_line_pipeline(device, &shader, &pipeline_layout, HDR_FORMAT);

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

        let (body_buf, body_cap) =
            alloc_instance_buf::<BodyInstance>(device, MIN_BUFFER_CAPACITY, "body");
        let (circle_buf, circle_cap) =
            alloc_instance_buf::<CircleInstance>(device, MIN_BUFFER_CAPACITY, "circle");
        let (line_buf, line_cap) =
            alloc_instance_buf::<LineInstance>(device, MIN_BUFFER_CAPACITY, "line");

        Self {
            body_pipeline,
            circle_pipeline,
            line_pipeline,
            screen_buf,
            lighting_buf,
            bind_group_screen,
            bind_group_layout_screen: screen_bgl,
            body_buf,
            body_cap,
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
        bodies: &[BodyInstance],
        circles: &[CircleInstance],
        lines: &[LineInstance],
    ) -> (u32, u32, u32) {
        queue.write_buffer(&self.screen_buf, 0, bytemuck::bytes_of(&screen));
        queue.write_buffer(&self.lighting_buf, 0, bytemuck::bytes_of(&lighting));

        if !bodies.is_empty() {
            ensure_instance_capacity::<BodyInstance>(
                device,
                queue,
                &mut self.body_buf,
                &mut self.body_cap,
                bodies.len() as u32,
                "body",
            );
            queue.write_buffer(&self.body_buf, 0, bytemuck::cast_slice(bodies));
        }

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

        (bodies.len() as u32, circles.len() as u32, lines.len() as u32)
    }

    fn draw(
        &self,
        pass: &mut wgpu::RenderPass<'_>,
        body_count: u32,
        circle_count: u32,
        line_count: u32,
    ) {
        // Layer order: lines (orbits) → annotation circles (apsides rings) →
        // bodies on top. Bodies are the visual focus and must not be
        // occluded by their own orbit markers.
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

        if body_count > 0 {
            pass.set_pipeline(&self.body_pipeline);
            pass.set_bind_group(0, &self.bind_group_screen, &[]);
            pass.set_vertex_buffer(0, self.body_buf.slice(..));
            pass.draw(0..6, 0..body_count);
        }
    }
}

// ── WgpuBackend ───────────────────────────────────────────────────────────────

pub struct WgpuBackend {
    bodies: Vec<BodyInstance>,
    circles: Vec<CircleInstance>,
    lines: Vec<LineInstance>,

    /// Global lighting state packed for the GPU. Rebuilt each frame via
    /// [`WgpuBackend::set_scene_lighting`] from a high-level
    /// [`SceneLighting`] so the shader only sees the GPU-ready shape.
    lighting: LightingUniform,

    gpu: Option<GpuResources>,
    trail: Option<TrailRenderer>,
    grid: Option<GridRenderer>,

    /// HDR scene colour target. Scene passes render into this linear-space
    /// texture; the tonemap pass composites it onto the swapchain.
    hdr: Option<HdrTarget>,
    tonemap: Option<TonemapPipeline>,

    /// GPU-side luminance metering for auto-exposure. Reduces the HDR
    /// target to a single `mean(L^p)` texel per frame; the CPU side
    /// reads it back through `exposure` below.
    luma_reducer: Option<LuminanceReducer>,
    /// EMA-smoothed exposure scale fed into the tonemap pipeline.
    pub exposure: ExposureState,
    /// Timestamp of the previous exposure tick — the EMA half-life
    /// needs wall-clock dt, not frame count, so adaptation feels the
    /// same at 30 fps and 240 fps.
    last_exposure_tick: Option<Instant>,

    pub trail_buffer: Option<Arc<apsis::core::trail::TrailBuffer>>,
    pub trail_visibility: Option<Vec<bool>>,
    pub center: [f32; 2],
    pub scale: f32,
    pub show_grid: bool,
    /// Visual style for trail rendering. Injected as a value object so
    /// presets swap atomically. See [`crate::render::TrailStylePreset`].
    pub trail_style: crate::render::TrailStyle,
}

impl WgpuBackend {
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            circles: Vec::new(),
            lines: Vec::new(),
            lighting: LightingUniform::default(),
            gpu: None,
            trail: None,
            grid: None,
            hdr: None,
            tonemap: None,
            luma_reducer: None,
            exposure: ExposureState::default(),
            last_exposure_tick: None,

            trail_buffer: None,
            trail_visibility: None,
            center: [0.0, 0.0],
            scale: 1.0,
            show_grid: true,
            trail_style: crate::render::TrailStylePreset::UniverseSandbox.style(1.5),
        }
    }

    pub fn begin(&mut self) {
        self.bodies.clear();
        self.circles.clear();
        self.lines.clear();
        // Reset lighting to the empty-scene default. The canvas layer
        // overwrites it later via `set_scene_lighting` if the frame has any
        // luminous bodies; otherwise the body shader falls back to pure
        // ambient × albedo + emissive.
        self.lighting = LightingUniform::default();
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

    /// Installs the per-frame scene lighting — sources, ambient floor,
    /// attenuation, and terminator-softening knob. The backend packs it
    /// into its GPU-uniform shape (sorting by intensity, clipping to
    /// [`crate::render::lighting::MAX_LIGHTS`], pre-squaring distances);
    /// callers pass the high-level [`SceneLighting`] and stay out of the
    /// byte layout.
    ///
    /// Call with `SceneLighting { lights: Vec::new(), .. }` on empty /
    /// dark systems — the shader's multi-light loop runs zero iterations
    /// and the body output collapses to `ambient_floor * albedo + emissive`,
    /// which is the correct fallback without a special-case flag.
    pub fn set_scene_lighting(&mut self, scene: SceneLighting) {
        self.lighting = LightingUniform::pack(&scene);
    }

    // ── Draw API ──────────────────────────────────────────────────────────────

    /// Submits a lit body — rendered via the sphere-from-disc pipeline with
    /// Lambertian diffuse from the primary light and unlit emissive on top.
    ///
    /// `screen_pos` is the body centre in logical pixels; `radius_screen` is
    /// its projected radius. `world_pos` carries the 3D world position used
    /// for the Lambert direction (`z = 0` for 2D top-down). `albedo` is the
    /// diffuse base colour; `emissive` the self-lit term (stars set this to
    /// their body colour and keep albedo near zero).
    pub fn draw_body(
        &mut self,
        screen_pos: [f32; 2],
        radius_screen: f32,
        world_pos: [f32; 3],
        albedo: [f32; 4],
        emissive: [f32; 4],
    ) {
        self.bodies.push(BodyInstance {
            center_screen: screen_pos,
            center_world: world_pos,
            radius_screen: radius_screen.max(0.5),
            albedo,
            emissive,
        });
    }

    /// Submits an annular ring (stroke) — flat, unlit. Used for annotations
    /// like apsides markers and selection rings.
    pub fn draw_circle_stroke(&mut self, pos: [f32; 2], radius: f32, width: f32, color: [u8; 4]) {
        let half = (width * 0.5).max(0.25);
        let outer = (radius + half).max(0.5);
        let inner = (radius - half).clamp(0.0, outer);

        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: outer,
            inner_radius: inner,
            color: rgba_u8_to_f32(color[0], color[1], color[2], color[3]),
        });
    }

    pub fn draw_line_segment(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: [u8; 4]) {
        // Sub-pixel line trick. Below 1 px wide the anti-alias band is
        // wider than the solid core, so the whole line smears into partial
        // alpha and reads as a faded ghost (worst in HDR+ACES, where low
        // mid-tones lose contrast against the black backdrop). Clamp the
        // rendered width to 1 px and scale alpha by the shortfall so total
        // coverage is preserved — a 0.8 px line becomes 1 px at 80% alpha,
        // visually identical but crisp instead of blurred.
        let alpha_scale = width.clamp(0.0, 1.0);
        let half_width = (width * 0.5).max(0.5);
        let a_scaled = (color[3] as f32 * alpha_scale).round().clamp(0.0, 255.0) as u8;

        self.lines.push(LineInstance {
            from,
            to,
            half_width,
            _pad: 0.0,
            color: rgba_u8_to_f32(color[0], color[1], color[2], a_scaled),
        });
    }

    // ── GPU setup ─────────────────────────────────────────────────────────────

    /// Creates the GPU pipelines on first use. Scene pipelines target the HDR
    /// offscreen format; only the tonemap pipeline depends on the swapchain
    /// format.
    pub fn ensure_gpu(&mut self, device: &wgpu::Device, swapchain_format: wgpu::TextureFormat) {
        if self.gpu.is_none() {
            self.gpu = Some(GpuResources::init(device));
        }
        if self.trail.is_none() {
            let gpu = self.gpu.as_ref().unwrap();
            self.trail =
                Some(TrailRenderer::new(device, &gpu.bind_group_layout_screen, HDR_FORMAT));
        }
        if self.grid.is_none() {
            self.grid = Some(GridRenderer::new(device, HDR_FORMAT));
        }
        if self.tonemap.is_none() {
            self.tonemap = Some(TonemapPipeline::new(device, swapchain_format));
        }
        if self.luma_reducer.is_none() {
            self.luma_reducer = Some(LuminanceReducer::new(device));
        }
    }

    // ── Render: split into scene (offscreen) + composite (tonemap) ───────────

    /// Records the scene into the HDR offscreen target. Runs in
    /// [`eframe::egui_wgpu::CallbackTrait::prepare`], where the encoder is
    /// available. The HDR texture is resized to match `physical_size` (canvas
    /// size in device pixels) before the pass begins.
    pub fn prepare_scene(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        encoder: &mut wgpu::CommandEncoder,
        screen: [f32; 2],
        viewport_min: [f32; 2],
        physical_size: [u32; 2],
        pixels_per_point: f32,
        swapchain_format: wgpu::TextureFormat,
        center: [f32; 2],
        scale: f32,
    ) {
        self.ensure_gpu(device, swapchain_format);
        self.center = center;
        self.scale = scale;

        // Keep the HDR target sized to the canvas.
        let hdr = self.hdr.get_or_insert_with(|| HdrTarget::new(device, physical_size));
        hdr.ensure_size(device, physical_size);

        // ── Auto-exposure tick (runs before tonemap upload) ──────────────────
        //
        // The reducer's reading describes the HDR frame submitted 1–2 frames
        // ago. We tick the EMA against that value and feed the resulting
        // scale into this frame's tonemap upload. The 1-frame lag is
        // imperceptible because the EMA already smooths over ~0.5 s of
        // adaptation; what matters is that the scale we upload *now* is
        // consistent with the reducer reading we observed *now*.
        if let Some(reducer) = self.luma_reducer.as_mut() {
            let now = Instant::now();
            let dt = self.last_exposure_tick.map(|t| (now - t).as_secs_f32()).unwrap_or(0.0);
            self.last_exposure_tick = Some(now);

            // `poll` returns the freshest reading that matured this frame;
            // if nothing matured, fall back to the last successful one so
            // the EMA keeps advancing toward the same target across the
            // empty frames. A completely cold start (no reading ever) just
            // leaves `current_scale` at its default.
            let raw = reducer.poll(device).or_else(|| reducer.last_reading());
            if let Some(mean_l_to_p) = raw {
                let soft_max = decode_reduced_texel(mean_l_to_p);
                self.exposure.tick(soft_max, dt);
            }
        }

        // Refresh the tonemap bind group if the HDR view was reallocated,
        // push the current exposure scale, and upload the uniform.
        if let Some(tm) = self.tonemap.as_mut() {
            tm.set_exposure(self.exposure.current_scale);
            tm.refresh_if_resized(device, hdr);
            tm.upload(queue);
        }

        let screen_uniform = ScreenUniform { size: screen, viewport_min };

        let (body_count, circle_count, line_count) = {
            let gpu = self.gpu.as_mut().unwrap();
            gpu.upload(
                device,
                queue,
                screen_uniform,
                self.lighting,
                &self.bodies,
                &self.circles,
                &self.lines,
            )
        };

        if self.show_grid {
            if let Some(grid) = &self.grid {
                // Grid fragment shader samples `frag_coord.xy` in the HDR
                // target's pixel space (canvas-local, physical pixels).
                // Camera state comes from the app in logical points inside the
                // full window, so translate by the canvas viewport and rescale
                // by pixels_per_point before upload. Without this the grid
                // lines and labels drift apart on the sub-pixel.
                let grid_center = [
                    (center[0] - viewport_min[0]) * pixels_per_point,
                    (center[1] - viewport_min[1]) * pixels_per_point,
                ];
                let grid_scale = scale * pixels_per_point;
                grid.upload(
                    queue,
                    grid_center,
                    grid_scale,
                    [physical_size[0] as f32, physical_size[1] as f32],
                );
            }
        }

        if let (Some(trail_renderer), Some(trail_buf)) =
            (self.trail.as_mut(), self.trail_buffer.as_deref())
        {
            trail_renderer.upload(
                device,
                queue,
                trail_buf,
                self.trail_visibility.as_deref(),
                center,
                scale,
                &self.trail_style,
            );
        }

        // ── Offscreen scene pass ─────────────────────────────────────────────
        // Wrapped in its own scope so the `RenderPass` is dropped before the
        // luminance reducer can record further passes on the same encoder.
        {
            let hdr_view = self.hdr.as_ref().unwrap().view();
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("scene::hdr_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: hdr_view,
                    resolve_target: None,
                    depth_slice: None,
                    ops: wgpu::Operations {
                        // Transparent clear: anywhere nothing is drawn, the
                        // tonemap composite lets the egui backdrop show
                        // through.
                        load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });

            // Scene layer order: grid → trails → primitives (circles + lines).
            if self.show_grid {
                if let Some(grid) = &self.grid {
                    grid.draw(&mut pass);
                }
            }

            if let Some(trail_renderer) = self.trail.as_ref() {
                if self.trail_buffer.is_some() {
                    let gpu = self.gpu.as_ref().unwrap();
                    trail_renderer.draw(&mut pass, &gpu.bind_group_screen);
                }
            }

            {
                let gpu = self.gpu.as_ref().unwrap();
                gpu.draw(&mut pass, body_count, circle_count, line_count);
            }
        }

        // ── Luminance reduce chain ───────────────────────────────────────────
        // Samples the freshly-written HDR target, reduces down to a single
        // `mean(L^p)` texel, and schedules an async readback. The result is
        // consumed in the *next* frame's exposure tick above.
        if self.exposure.enabled {
            if let Some(reducer) = self.luma_reducer.as_mut() {
                let hdr = self.hdr.as_ref().unwrap();
                reducer.encode(device, queue, encoder, hdr);
            }
        }
    }

    /// Composites the HDR scene target onto the supplied swapchain pass via
    /// the tonemap pipeline. Runs in
    /// [`eframe::egui_wgpu::CallbackTrait::paint`].
    pub fn composite(&self, pass: &mut wgpu::RenderPass<'_>) {
        if let Some(tm) = self.tonemap.as_ref() {
            tm.draw(pass);
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

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
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
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
    // center(2) + outer_r(1) + inner_r(1) + color(4)
    let attrs = wgpu::vertex_attr_array![
        0 => Float32x2,   // center
        1 => Float32,     // outer_radius
        2 => Float32,     // inner_radius
        3 => Float32x4    // color
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

fn build_body_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    layout: &wgpu::PipelineLayout,
    format: wgpu::TextureFormat,
) -> wgpu::RenderPipeline {
    // center_screen(2) + center_world(3) + radius_screen(1) + albedo(4) + emissive(4).
    // Tight packing: offsets 0, 8, 20, 24, 40 — matches #[repr(C)] BodyInstance.
    let attrs = wgpu::vertex_attr_array![
        0 => Float32x2,   // center_screen
        1 => Float32x3,   // center_world
        2 => Float32,     // radius_screen
        3 => Float32x4,   // albedo
        4 => Float32x4    // emissive
    ];

    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("grav_sim::body_pipeline"),
        layout: Some(layout),
        vertex: wgpu::VertexState {
            module: shader,
            entry_point: Some("vs_body"),
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: size_of::<BodyInstance>() as u64,
                step_mode: wgpu::VertexStepMode::Instance,
                attributes: &attrs,
            }],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: shader,
            entry_point: Some("fs_body"),
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

/// Multi-light scene configuration. Matches [`crate::render::lighting::LightingUniform`]
/// byte-for-byte — see that type for field semantics.
///
/// 3D-ready: every light position is `vec3<f32>`; 2D callers set `z = 0`,
/// a future 3D camera populates the full vector with no shader change.
/// The array size `4` mirrors `MAX_LIGHTS`; bumping one side without the
/// other will be caught at pipeline creation.
struct PackedLight {
    world_pos: vec3<f32>,
    intensity: f32,
};

struct LightingUniform {
    lights:        array<PackedLight, 4>,
    num_lights:    u32,
    ambient_floor: f32,
    r_ref_sq:      f32,
    bias_sq:       f32,
    wrap:          f32,
    _pad0:         f32,
    _pad1:         f32,
    _pad2:         f32,
};

@group(0) @binding(0) var<uniform> screen:   ScreenUniform;
@group(0) @binding(1) var<uniform> lighting: LightingUniform;

fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let local = p - screen.viewport_min;
    let x =  (local.x / screen.size.x) * 2.0 - 1.0;
    let y = -(local.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

// ── CIRCLES (flat annotations) ───────────────────────────────────────────────
// Used for orbit apsides markers, selection rings, etc. — never lit.

struct CircleVarying {
    @builtin(position) clip_pos:    vec4<f32>,
    @location(0)       local:       vec2<f32>,
    @location(1)       inner_ratio: f32,
    @location(2)       color:       vec4<f32>,
};

@vertex
fn vs_circle(
    @builtin(vertex_index) vi:           u32,
    @location(0)           center:       vec2<f32>,
    @location(1)           outer_radius: f32,
    @location(2)           inner_radius: f32,
    @location(3)           color:        vec4<f32>,
) -> CircleVarying {
    var quad = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );

    let local     = quad[vi];
    let world_pos = center + local * outer_radius;

    var out: CircleVarying;
    out.clip_pos    = to_ndc(world_pos);
    out.local       = local;
    out.inner_ratio = select(0.0, inner_radius / outer_radius, outer_radius > 0.0);
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

    let alpha = in.color.a * outer * inner_mask;
    if alpha <= 0.001 { discard; }
    return vec4<f32>(in.color.rgb, alpha);
}

// ── BODIES (sphere-from-disc + multi-light + emissive) ───────────────────────
// Bodies expand into a screen-aligned quad per instance. The fragment shader
// reconstructs a unit-sphere surface normal from the quad's (u, v) local
// coordinate — world X/Y align with screen X/Y, Z is out of the screen toward
// the camera — and loops over every registered light, summing wrap-softened
// Lambert contributions attenuated by r_ref²/(r²+bias²). Emissive is added
// unattenuated on top so stars self-illuminate; a multiplicative ambient
// floor keeps the dark side from vanishing into the HDR black.

struct BodyVarying {
    @builtin(position) clip_pos:      vec4<f32>,
    /// Quad-local coordinate in [-1, 1]². Doubles as the XY sphere normal.
    @location(0)       local:         vec2<f32>,
    /// Body center in world coordinates. Fragments compute their own
    /// world position as `center_world + vec3(local, 0) * scale` — but
    /// since the body radius is negligible versus inter-body distance,
    /// we treat the whole disc as sitting at the center for lighting.
    /// That keeps the per-fragment cost to a scalar divide per light.
    @location(1)       center_world:  vec3<f32>,
    @location(2)       albedo:        vec4<f32>,
    @location(3)       emissive:      vec4<f32>,
    /// Screen radius in logical pixels — used in the fragment to blend
    /// toward flat shading on sub-pixel bodies (see `fs_body`).
    @location(4)       radius_screen: f32,
};

@vertex
fn vs_body(
    @builtin(vertex_index) vi:             u32,
    @location(0)           center_screen:  vec2<f32>,
    @location(1)           center_world:   vec3<f32>,
    @location(2)           radius_screen:  f32,
    @location(3)           albedo:         vec4<f32>,
    @location(4)           emissive:       vec4<f32>,
) -> BodyVarying {
    var quad = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );

    let local    = quad[vi];
    let screen_p = center_screen + local * radius_screen;

    var out: BodyVarying;
    out.clip_pos      = to_ndc(screen_p);
    out.local         = local;
    out.center_world  = center_world;
    out.albedo        = albedo;
    out.emissive      = emissive;
    out.radius_screen = radius_screen;
    return out;
}

@fragment
fn fs_body(in: BodyVarying) -> @location(0) vec4<f32> {
    let r2 = dot(in.local, in.local);
    let r  = sqrt(r2);
    let aa = fwidth(r);
    let mask = 1.0 - smoothstep(1.0 - aa, 1.0 + aa, r);
    if mask <= 0.001 { discard; }

    // Sphere-from-disc normal. The camera looks down -Z in world, so the
    // visible hemisphere has n.z > 0. `max(0, …)` guards the corner near
    // r = 1 where floating-point noise could drive the argument slightly
    // negative and turn `sqrt` into NaN.
    let nz = sqrt(max(0.0, 1.0 - r2));
    let n  = vec3<f32>(in.local, nz);

    // Accumulate per-light contributions. The loop over num_lights runs
    // zero iterations on empty / dark systems, so the body collapses to
    // ambient × albedo + emissive without needing a has_light flag.
    var diffuse_total = 0.0;
    for (var i: u32 = 0u; i < lighting.num_lights; i = i + 1u) {
        let L       = lighting.lights[i];
        let to_l    = L.world_pos - in.center_world;
        let d2      = dot(to_l, to_l);
        // Pre-squared reference distance + bias: attenuation = r_ref²/(r²+bias²).
        // The bias prevents the singularity when a body overlaps its own light
        // (e.g. a star shading itself) and gives an artistic knob for soft falloff.
        let att     = lighting.r_ref_sq / (d2 + lighting.bias_sq);
        let inv_len = inverseSqrt(max(d2, 1e-20));
        let dir     = to_l * inv_len;

        let ndotl   = dot(n, dir);
        let lambert = max(ndotl, 0.0);
        // Half-Lambert wrap: ((n·L + 1) / 2)². Softens the terminator so
        // back-facing fragments still pick up a little light — physically
        // motivated for dusty / atmospheric bodies, aesthetically for
        // avoiding hard black crescents that eat the silhouette.
        let wrapped = ((ndotl + 1.0) * 0.5) * ((ndotl + 1.0) * 0.5);
        let shading = mix(lambert, wrapped, lighting.wrap);

        diffuse_total = diffuse_total + shading * att * L.intensity;
    }

    // Multiplicative ambient floor: mix(floor, 1, diffuse). Small values
    // (~0.05) keep the unlit side readable in HDR black without washing
    // out the lit side — additive ambient would over-brighten the day
    // hemisphere since it's already at 1.0.
    let lit_factor = mix(lighting.ambient_floor, 1.0, clamp(diffuse_total, 0.0, 1.0));

    // Sub-pixel resolution fallback.
    // When a body is only a few pixels wide, the Lambert terminator sweeps
    // across most of the disc and the unlit hemisphere reads as "half-moon
    // cutout" or "missing body" against the black backdrop — the eye can't
    // resolve the phase at that scale anyway. Blend `lit_factor` toward 1.0
    // (flat-shaded full disc) as the screen radius shrinks below ~6 px, so
    // distant planets still register as solid coloured dots while the
    // beautifully-phased large discs keep their crescents intact.
    let flat_weight = 1.0 - smoothstep(2.5, 6.0, in.radius_screen);
    let effective   = mix(lit_factor, 1.0, flat_weight);

    let rgb   = in.albedo.rgb * effective + in.emissive.rgb;
    let alpha = in.albedo.a * mask;
    return vec4<f32>(rgb, alpha);
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

    // Overlap neighbours by ~0.5 logical px on each end. A polyline is a
    // sequence of butt-capped rectangles; where two segments meet at an
    // angle, the butt ends leave a small triangular gap that reads as a
    // dashed line on thin strokes. Extending the quad slightly past each
    // endpoint makes adjacent segments overlap by ~1 px, closing the gap.
    // Capped at len/2 so segments below 1 px long don't flip.
    let overlap = min(0.5, len * 0.49);
    let x_along = corner.x * (len + overlap * 2.0) - overlap;
    let pos     = p0 + tangent * x_along + normal * (half_width * corner.y);

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
