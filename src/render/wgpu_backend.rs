//! WGPU primitive renderer for the gravity simulation canvas.
//!
//! # Overview
//!
//! [`WgpuBackend`] collects draw commands during a frame and submits them as a
//! single [`egui_wgpu::Callback`] on [`RenderBackend::end`].  Circles and line
//! segments are rendered via GPU instancing — one draw call per primitive type
//! per frame.
//!
//! # Buffer strategy
//!
//! All persistent GPU state lives in [`GpuResources`], which is stored in
//! [`egui_wgpu::CallbackResources`] and survives across frames.  Per-frame
//! instance data is uploaded with [`wgpu::Queue::write_buffer`] — **no GPU
//! allocation occurs on the hot path** once the buffer is large enough.
//! Buffers grow on demand using power-of-two doubling (see
//! [`ensure_instance_capacity`]).
//!
//! # Prepare → paint data flow
//!
//! [`WgpuPrimitivesCallback`] implements [`egui_wgpu::CallbackTrait`].  Both
//! `prepare` and `paint` receive `&self`, so the instance counts produced in
//! `prepare` are forwarded to `paint` through [`AtomicU32`] fields — lock-free
//! and zero-overhead for the single-writer / single-reader pattern here.

use std::mem::size_of;
use std::sync::atomic::{AtomicU32, Ordering};

use bytemuck::{Pod, Zeroable};
use eframe::egui::{self, Rect};
use eframe::egui_wgpu::{self, CallbackTrait, ScreenDescriptor};

use crate::render::RenderBackend;

// ── Constants ─────────────────────────────────────────────────────────────────

/// Minimum (and initial) instance-buffer capacity, in element count.
///
/// Buffers are never created smaller than this.  Growth beyond this point uses
/// power-of-two doubling to amortise GPU-allocation cost.
const MIN_BUFFER_CAPACITY: u32 = 256;

// ── GPU-mapped data types ─────────────────────────────────────────────────────

/// Screen dimensions sent to the vertex shader once per frame.
///
/// Values are in *logical* pixels (`physical_pixels / pixels_per_point`) so
/// the NDC mapping is DPI-independent and consistent with egui's layout space.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ScreenUniform {
    /// `[logical_width, logical_height]`.
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

/// Per-instance data for a filled or stroked disc.
///
/// `inner_radius = 0.0` produces a solid disc; a positive value produces an
/// annular ring (stroke) with the given inner cutout.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CircleInstance {
    /// Centre in logical screen pixels.
    center: [f32; 2],
    /// Outer radius in logical pixels.  Clamped to `≥ 0.5` by the backend.
    outer_radius: f32,
    /// Inner radius in logical pixels.  `0.0` = filled disc.
    inner_radius: f32,
    _pad: f32,
    /// RGBA colour, each component in `[0.0, 1.0]`.
    color: [f32; 4],
}

/// Per-instance data for a line segment rendered as a screen-aligned quad.
///
/// The vertex shader extrudes the segment by `half_width` on each side of the
/// line axis to produce a thick, filled rectangle.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct LineInstance {
    /// Segment start in logical screen pixels.
    from: [f32; 2],
    /// Segment end in logical screen pixels.
    to: [f32; 2],
    /// Half the stroke width in logical pixels.  Clamped to `≥ 0.25`.
    half_width: f32,
    _pad: f32,
    /// RGBA colour, each component in `[0.0, 1.0]`.
    color: [f32; 4],
}

// ── Persistent GPU resources ──────────────────────────────────────────────────

/// All GPU-side state that persists across frames.
///
/// Stored in [`egui_wgpu::CallbackResources`] and initialised lazily by
/// [`WgpuPrimitivesCallback`] on the first call to `prepare`.
///
/// ## Invariants
///
/// - `screen_buf` is **never reallocated**; the bind group therefore remains
///   valid for the entire application lifetime.
/// - `circle_buf` / `line_buf` may be **replaced** (full realloc, old handle
///   dropped) when [`ensure_instance_capacity`] triggers growth.  The bind
///   group is unaffected because it references only `screen_buf`.
/// - Buffer capacities are always powers of two ≥ [`MIN_BUFFER_CAPACITY`].
struct GpuResources {
    // ── Render pipelines ─────────────────────────────────────────────── //
    circle_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,

    // ── Screen uniform ────────────────────────────────────────────────── //
    /// 16-byte uniform; overwritten every frame via `write_buffer`.
    screen_buf: wgpu::Buffer,
    /// Permanent bind group pointing at `screen_buf`.
    bind_group: wgpu::BindGroup,

    // ── Instance buffers ──────────────────────────────────────────────── //
    /// Vertex buffer for [`CircleInstance`] data.
    circle_buf: wgpu::Buffer,
    /// Allocated slot count in `circle_buf`.
    circle_cap: u32,

    /// Vertex buffer for [`LineInstance`] data.
    line_buf: wgpu::Buffer,
    /// Allocated slot count in `line_buf`.
    line_cap: u32,
}

impl GpuResources {
    /// Creates all GPU resources.
    ///
    /// Called once; the result is inserted into
    /// [`egui_wgpu::CallbackResources`] and reused for every subsequent frame.
    fn init(device: &wgpu::Device) -> Self {
        // All intermediate handles (bgl, shader, pipeline_layout) are ref-counted
        // on the GPU side; dropping them here after pipeline creation is safe.
        let bgl = build_bind_group_layout(device);
        let shader = build_shader(device);
        let pipeline_layout = build_pipeline_layout(device, &bgl);

        let circle_pipeline = build_circle_pipeline(device, &shader, &pipeline_layout);
        let line_pipeline = build_line_pipeline(device, &shader, &pipeline_layout);

        let screen_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("grav_sim::screen_uniform"),
            size: size_of::<ScreenUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // The bind group references `screen_buf` by object identity, so
        // subsequent `write_buffer` calls do not invalidate it.
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("grav_sim::bind_group"),
            layout: &bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: screen_buf.as_entire_binding(),
            }],
        });

        let (circle_buf, circle_cap) =
            alloc_instance_buf::<CircleInstance>(device, MIN_BUFFER_CAPACITY, "grav_sim::circle");
        let (line_buf, line_cap) =
            alloc_instance_buf::<LineInstance>(device, MIN_BUFFER_CAPACITY, "grav_sim::line");

        Self {
            circle_pipeline,
            line_pipeline,
            screen_buf,
            bind_group,
            circle_buf,
            circle_cap,
            line_buf,
            line_cap,
        }
    }

    /// Uploads one frame's instance data and the screen dimensions.
    ///
    /// Writes `screen` to the persistent uniform buffer (16 bytes — negligible
    /// bandwidth) and copies instance slices to the vertex buffers, growing
    /// them with power-of-two doubling when the current capacity is exceeded.
    ///
    /// Returns `(circle_count, line_count)` for use by [`draw`](Self::draw).
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
                &mut self.circle_buf,
                &mut self.circle_cap,
                circles.len() as u32,
                "grav_sim::circle",
            );
            queue.write_buffer(&self.circle_buf, 0, bytemuck::cast_slice(circles));
        }

        if !lines.is_empty() {
            ensure_instance_capacity::<LineInstance>(
                device,
                &mut self.line_buf,
                &mut self.line_cap,
                lines.len() as u32,
                "grav_sim::line",
            );
            queue.write_buffer(&self.line_buf, 0, bytemuck::cast_slice(lines));
        }

        (circles.len() as u32, lines.len() as u32)
    }

    /// Issues GPU draw calls for the given instance counts.
    ///
    /// Lines are submitted before circles so that body discs render on top of
    /// trail segments when both occupy the same screen area.
    fn draw<'rp>(
        &'rp self,
        pass: &mut wgpu::RenderPass<'rp>,
        circle_count: u32,
        line_count: u32,
    ) {
        if line_count > 0 {
            pass.set_pipeline(&self.line_pipeline);
            pass.set_bind_group(0, &self.bind_group, &[]);
            pass.set_vertex_buffer(0, self.line_buf.slice(..));
            pass.draw(0..6, 0..line_count);
        }

        if circle_count > 0 {
            pass.set_pipeline(&self.circle_pipeline);
            pass.set_bind_group(0, &self.bind_group, &[]);
            pass.set_vertex_buffer(0, self.circle_buf.slice(..));
            pass.draw(0..6, 0..circle_count);
        }
    }
}

// ── Per-frame paint callback ──────────────────────────────────────────────────

/// One-frame paint callback registered with egui's WGPU integration.
///
/// Carries the CPU-side instance slices for the current frame.  GPU resources
/// are **not** recreated here — they live in [`GpuResources`].
///
/// ## Prepare → paint handoff
///
/// [`egui_wgpu::CallbackTrait`] gives both `prepare` and `paint` only a shared
/// reference (`&self`).  The instance counts produced in `prepare` are passed
/// to `paint` via [`AtomicU32`] fields with `Relaxed` ordering, which is
/// correct because egui_wgpu guarantees `prepare` completes before `paint`.
struct WgpuPrimitivesCallback {
    circles: Vec<CircleInstance>,
    lines: Vec<LineInstance>,
    /// Written by `prepare`, read by `paint`.
    circle_count: AtomicU32,
    /// Written by `prepare`, read by `paint`.
    line_count: AtomicU32,
}

impl CallbackTrait for WgpuPrimitivesCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        screen_descriptor: &ScreenDescriptor,
        _encoder: &mut wgpu::CommandEncoder,
        resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        if resources.get::<GpuResources>().is_none() {
            resources.insert(GpuResources::init(device));
        }
        let gpu = resources
            .get_mut::<GpuResources>()
            .expect("GpuResources should be present after lazy init");

        let screen = ScreenUniform::from_descriptor(screen_descriptor);
        let (cc, lc) = gpu.upload(device, queue, screen, &self.circles, &self.lines);

        self.circle_count.store(cc, Ordering::Relaxed);
        self.line_count.store(lc, Ordering::Relaxed);

        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        pass: &mut wgpu::RenderPass<'static>,
        resources: &egui_wgpu::CallbackResources,
    ) {
        let Some(gpu) = resources.get::<GpuResources>() else {
            return;
        };

        // SAFETY: `GpuResources` is stored inside `CallbackResources`, which
        // enforces a `'static` bound on all inserted types and is itself
        // application-scoped.  The reference lifetime is artificially narrowed
        // to that of `resources`; the underlying allocation lives for the
        // entire application lifetime.  We need `'static` here so that buffer
        // slices derived from `gpu` satisfy `RenderPass<'static>::set_vertex_buffer`.
        let gpu: &'static GpuResources = unsafe { &*(gpu as *const GpuResources) };

        gpu.draw(
            pass,
            self.circle_count.load(Ordering::Relaxed),
            self.line_count.load(Ordering::Relaxed),
        );
    }
}

// ── Public interface ──────────────────────────────────────────────────────────

/// WGPU-backed [`RenderBackend`] for the simulation canvas.
///
/// Accumulates draw commands into CPU-side `Vec`s and registers a single
/// [`egui_wgpu::Callback`] on [`end`](Self::end).  All GPU work is deferred to
/// the callback's `prepare` / `paint` cycle.
///
/// # Example
///
/// ```ignore
/// let mut backend = WgpuBackend::new(ui, canvas_rect);
/// backend.begin();
/// backend.draw_circle([cx, cy], radius, [r, g, b]);
/// backend.draw_line_segment([x0, y0], [x1, y1], width, [r, g, b, a]);
/// backend.end();
/// ```
pub struct WgpuBackend<'a> {
    ui: &'a egui::Ui,
    rect: Rect,
    circles: Vec<CircleInstance>,
    lines: Vec<LineInstance>,
}

impl<'a> WgpuBackend<'a> {
    /// Creates a new backend targeting `rect` within `ui`.
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
    /// Clears accumulated draw commands in preparation for a new frame.
    fn begin(&mut self) {
        self.circles.clear();
        self.lines.clear();
    }

    /// Queues a solid filled disc.
    fn draw_circle(&mut self, pos: [f32; 2], radius: f32, color: [u8; 3]) {
        self.circles.push(CircleInstance {
            center: pos,
            outer_radius: radius.max(0.5),
            inner_radius: 0.0,
            _pad: 0.0,
            color: rgba_u8_to_f32(color[0], color[1], color[2], 255),
        });
    }

    /// Queues an annular ring (stroke around a disc).
    fn draw_circle_stroke(&mut self, pos: [f32; 2], radius: f32, width: f32, color: [u8; 4]) {
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

    /// Queues a screen-aligned thick line segment.
    fn draw_line_segment(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: [u8; 4]) {
        self.lines.push(LineInstance {
            from,
            to,
            half_width: (width * 0.5).max(0.25),
            _pad: 0.0,
            color: rgba_u8_to_f32(color[0], color[1], color[2], color[3]),
        });
    }

    /// Flushes all queued commands as a single egui WGPU paint callback.
    ///
    /// The accumulated instance vecs are *moved into* the callback (zero-copy);
    /// this backend's buffers are left empty after the call.
    fn end(&mut self) {
        let callback = egui_wgpu::Callback::new_paint_callback(
            self.rect,
            WgpuPrimitivesCallback {
                circles: std::mem::take(&mut self.circles),
                lines: std::mem::take(&mut self.lines),
                circle_count: AtomicU32::new(0),
                line_count: AtomicU32::new(0),
            },
        );
        self.ui.painter().add(callback);
    }
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Converts a `[u8; 4]` colour to `[f32; 4]` with a simple `/ 255.0` scale.
///
/// No gamma linearisation is applied; this matches egui's convention for the
/// `Bgra8Unorm` render target used by the simulation canvas.
#[inline]
fn rgba_u8_to_f32(r: u8, g: u8, b: u8, a: u8) -> [f32; 4] {
    [
        r as f32 / 255.0,
        g as f32 / 255.0,
        b as f32 / 255.0,
        a as f32 / 255.0,
    ]
}

/// Allocates a new vertex buffer sized for `capacity` elements of type `T`.
///
/// The buffer is created with `VERTEX | COPY_DST` so it can be updated every
/// frame via [`wgpu::Queue::write_buffer`] without reallocation.
fn alloc_instance_buf<T: Pod>(
    device: &wgpu::Device,
    capacity: u32,
    label: &str,
) -> (wgpu::Buffer, u32) {
    let buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some(label),
        size: capacity as u64 * size_of::<T>() as u64,
        usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    (buf, capacity)
}

/// Ensures `buf` can hold at least `needed` elements of type `T`.
///
/// If the current capacity `cap` is insufficient, the buffer is reallocated
/// with `max(needed.next_power_of_two(), MIN_BUFFER_CAPACITY)` slots.  The
/// old buffer handle is dropped immediately.  No data is preserved or copied.
fn ensure_instance_capacity<T: Pod>(
    device: &wgpu::Device,
    buf: &mut wgpu::Buffer,
    cap: &mut u32,
    needed: u32,
    label: &str,
) {
    if needed <= *cap {
        return;
    }
    let new_cap = needed.next_power_of_two().max(MIN_BUFFER_CAPACITY);
    (*buf, *cap) = alloc_instance_buf::<T>(device, new_cap, label);
}

// ── Pipeline / shader builders (called once during GpuResources::init) ────────

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
) -> wgpu::RenderPipeline {
    // centre(2) + outer_r(1) + inner_r(1) + pad(1) + color(4) = 9 floats
    let attrs = wgpu::vertex_attr_array![
        0 => Float32x2, // center
        1 => Float32,   // outer_radius
        2 => Float32,   // inner_radius
        3 => Float32,   // _pad
        4 => Float32x4  // color
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

fn build_line_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    layout: &wgpu::PipelineLayout,
) -> wgpu::RenderPipeline {
    // from(2) + to(2) + half_width(1) + pad(1) + color(4) = 10 floats
    let attrs = wgpu::vertex_attr_array![
        0 => Float32x2, // from
        1 => Float32x2, // to
        2 => Float32,   // half_width
        3 => Float32,   // _pad
        4 => Float32x4  // color
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
                format: wgpu::TextureFormat::Bgra8Unorm,
                blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                write_mask: wgpu::ColorWrites::ALL,
            })],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        }),
        // The vertex shader generates a 2-triangle quad (6 vertices).
        // TriangleList is required to rasterise it as a filled rectangle.
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

// ── Shared utility ────────────────────────────────────────────────────────── //

/// Maps a logical-pixel position to WGPU clip space (NDC, Y-up).
fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let x =  (p.x / screen.size.x) * 2.0 - 1.0;
    let y = -(p.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

// ── Circle pipeline ───────────────────────────────────────────────────────── //

struct CircleVarying {
    @builtin(position) clip_pos:   vec4<f32>,
    @location(0)       local:      vec2<f32>,
    @location(1)       inner_ratio: f32,
    @location(2)       color:      vec4<f32>,
};

@vertex
fn vs_circle(
    @builtin(vertex_index) vi:    u32,
    @location(0) center:          vec2<f32>,
    @location(1) outer_radius:    f32,
    @location(2) inner_radius:    f32,
    @location(3) _pad:            f32,
    @location(4) color:           vec4<f32>,
) -> CircleVarying {
    // Unit quad in [-1, 1]², two triangles, CCW winding.
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
    let r2 = dot(in.local, in.local);
    // Discard outside the outer circle or inside the inner cutout.
    if r2 > 1.0 || r2 < in.inner_ratio * in.inner_ratio {
        discard;
    }
    return in.color;
}

// ── Line pipeline ─────────────────────────────────────────────────────────── //

struct LineVarying {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0)       color:    vec4<f32>,
};

@vertex
fn vs_line(
    @builtin(vertex_index) vi: u32,
    @location(0) p0:           vec2<f32>,
    @location(1) p1:           vec2<f32>,
    @location(2) half_width:   f32,
    @location(3) _pad:         f32,
    @location(4) color:        vec4<f32>,
) -> LineVarying {
    let dir     = p1 - p0;
    let len     = max(length(dir), 1e-5);
    let tangent = dir / len;
    let normal  = vec2<f32>(-tangent.y, tangent.x);

    // Quad corners: u in [0,1] along axis, v in [-1,1] across axis.
    // Six vertices = two CCW triangles.
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
    out.color    = color;
    return out;
}

@fragment
fn fs_line(in: LineVarying) -> @location(0) vec4<f32> {
    return in.color;
}
"#;
