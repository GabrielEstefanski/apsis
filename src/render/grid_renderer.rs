use bytemuck::{Pod, Zeroable};
use std::mem::size_of;

// ── Uniform ───────────────────────────────────────────────────────────────────

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct GridUniform {
    screen_size: [f32; 2],
    center: [f32; 2],
    scale: f32,
    _pad: [f32; 3],
}

// ── GridRenderer ──────────────────────────────────────────────────────────────

pub struct GridRenderer {
    pipeline: wgpu::RenderPipeline,
    uniform_buf: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
}

impl GridRenderer {
    pub fn new(device: &wgpu::Device, format: wgpu::TextureFormat) -> Self {
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("grid::bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("grid::shader"),
            source: wgpu::ShaderSource::Wgsl(GRID_SHADER.into()),
        });

        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("grid::layout"),
            bind_group_layouts: &[Some(&bgl)],
            immediate_size: 0,
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("grid::pipeline"),
            layout: Some(&layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_grid"),
                buffers: &[],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_grid"),
                targets: &[Some(wgpu::ColorTargetState {
                    format,
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

        let uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("grid::uniform"),
            size: size_of::<GridUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("grid::bg"),
            layout: &bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buf.as_entire_binding(),
            }],
        });

        Self {
            pipeline,
            uniform_buf,
            bind_group,
        }
    }

    pub fn upload(&self, queue: &wgpu::Queue, center: [f32; 2], scale: f32, screen: [f32; 2]) {
        let u = GridUniform {
            screen_size: screen,
            center,
            scale,
            _pad: [0.0; 3],
        };
        queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&u));
    }

    pub fn draw(&self, pass: &mut wgpu::RenderPass<'_>) {
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.draw(0..6, 0..1);
    }

    /// Overlay coordinate labels on top of the GPU grid using egui.
    ///
    /// Call this every frame after [`draw`], inside your egui paint callback
    /// or immediate-mode panel that covers the simulation viewport.
    ///
    /// Labels are drawn at the "medium" grid level (every 5 fine lines) so
    /// they are sparse enough not to clutter the view. The label for the
    /// origin is suppressed to avoid overlap with the axis intersection.
    ///
    /// `rect` must be the egui [`Rect`] that matches the wgpu viewport so
    /// that pixel coordinates are consistent between the two renderers.
    #[cfg(feature = "egui")]
    pub fn draw_labels(
        &self,
        ui: &egui::Ui,
        center: [f32; 2],
        scale: f32,
        screen: [f32; 2],
        rect: egui::Rect,
    ) {
        let painter = ui.painter_at(rect);
        let font = egui::FontId::monospace(10.0);
        let color = egui::Color32::from_rgba_unmultiplied(180, 180, 220, 180);

        let step = nice_number(60.0 / scale) * 5.0;

        // ── X-axis labels (placed just below the x-axis, clamped to viewport) ──
        let y_label = center[1].clamp(rect.min.y + 16.0, rect.max.y - 16.0);
        let x_world_min = (rect.min.x - center[0]) / scale;
        let x_world_max = (rect.max.x - center[0]) / scale;

        let ix_min = (x_world_min / step).floor() as i32;
        let ix_max = (x_world_max / step).ceil() as i32;

        for ix in ix_min..=ix_max {
            let wx = ix as f32 * step;
            if wx.abs() < step * 0.1 {
                continue;
            } // suppress origin
            let sx = center[0] + wx * scale;
            painter.text(
                egui::pos2(sx, y_label + 4.0),
                egui::Align2::CENTER_TOP,
                format_coord(wx),
                font.clone(),
                color,
            );
        }

        // ── Y-axis labels (placed just left of the y-axis, clamped to viewport) ──
        let x_label = center[0].clamp(rect.min.x + 36.0, rect.max.x - 8.0);
        let y_world_min = (rect.min.y - center[1]) / scale;
        let y_world_max = (rect.max.y - center[1]) / scale;

        let iy_min = (y_world_min / step).floor() as i32;
        let iy_max = (y_world_max / step).ceil() as i32;

        for iy in iy_min..=iy_max {
            let wy = iy as f32 * step;
            if wy.abs() < step * 0.1 {
                continue;
            } // suppress origin
            let sy = center[1] + wy * scale;
            painter.text(
                egui::pos2(x_label - 4.0, sy),
                egui::Align2::RIGHT_CENTER,
                format_coord(wy),
                font.clone(),
                color,
            );
        }

        // ── Origin label ──────────────────────────────────────────────────────
        painter.text(
            egui::pos2(
                center[0].clamp(rect.min.x + 20.0, rect.max.x - 8.0) - 4.0,
                center[1].clamp(rect.min.y + 8.0, rect.max.y - 8.0) + 4.0,
            ),
            egui::Align2::RIGHT_TOP,
            "0",
            font,
            color,
        );
    }
}

// ── Coordinate formatting ─────────────────────────────────────────────────────

/// Format a world coordinate for display as a grid label.
///
/// Uses compact scientific notation (`1.2e4`) for very large or very small
/// values, plain integers when the value is whole, and two decimal places
/// otherwise. This keeps labels readable across many orders of magnitude,
/// which is common in astronomical simulations.
fn format_coord(v: f32) -> String {
    let abs = v.abs();
    if abs == 0.0 {
        return "0".into();
    }
    if abs >= 1.0e4 || abs < 1.0e-2 {
        let exp = abs.log10().floor() as i32;
        let mantissa = v / 10_f32.powi(exp);
        return if (mantissa - mantissa.round()).abs() < 0.05 {
            format!("{}e{exp}", mantissa.round() as i32)
        } else {
            format!("{mantissa:.1}e{exp}")
        };
    }
    if (v - v.round()).abs() < 1.0e-3 {
        return format!("{}", v.round() as i32);
    }
    format!("{v:.2}")
}

/// Round `raw` up to the nearest "nice" number: 1, 2, 5, 10, 20, 50, …
///
/// This is the standard Wilkinson / Heckbert algorithm used by most
/// scientific plotting libraries to choose axis tick spacing.
fn nice_number(raw: f32) -> f32 {
    let e = raw.abs().log10().floor();
    let base = 10_f32.powf(e);
    let frac = raw / base;
    if frac < 1.5 {
        base
    } else if frac < 3.5 {
        base * 2.0
    } else if frac < 7.5 {
        base * 5.0
    } else {
        base * 10.0
    }
}

// ── WGSL Shader ───────────────────────────────────────────────────────────────

const GRID_SHADER: &str = r#"

struct GridUniform {
    screen_size : vec2<f32>,
    center      : vec2<f32>,
    scale       : f32,
    _pad0       : f32,
    _pad1       : f32,
    _pad2       : f32,
};

@group(0) @binding(0)
var<uniform> u : GridUniform;

// ── Vertex: procedural full-screen quad, no vertex buffer ─────────────────── //

@vertex
fn vs_grid(@builtin(vertex_index) vi : u32) -> @builtin(position) vec4<f32> {
    var verts = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );
    return vec4<f32>(verts[vi], 0.0, 1.0);
}

// ── Helpers ───────────────────────────────────────────────────────────────── //

/// Antialiased coverage of a grid line at coordinate `v` with spacing `s`.
/// Uses the screen-space derivative (fwidth) of the normalised coordinate so
/// line width is always exactly one pixel regardless of zoom level.
/// Returns 0..1 where 1 is the centre of the line.
fn grid_aa(v : f32, s : f32) -> f32 {
    let c    = v / s;
    let d    = fwidth(c);
    let dist = abs(fract(c + 0.5) - 0.5);
    return 1.0 - smoothstep(0.0, d, dist);
}

/// Maximum coverage across both grid axes for a 2-D grid.
fn grid_line(world : vec2<f32>, s : f32) -> f32 {
    return max(grid_aa(world.x, s), grid_aa(world.y, s));
}

/// Per-axis color and coverage for the principal axes (x = 0, y = 0).
///
/// Follows the scientific / CAD convention:
///   X-axis (y = 0) → muted red
///   Y-axis (x = 0) → muted green
///
/// Returns vec4(rgb, alpha). Alpha is 0 when neither axis is close.
fn axis_color(world : vec2<f32>) -> vec4<f32> {
    let px        = abs(world) * u.scale;
    let threshold = 1.2;
    // ax: coverage of the Y-axis (vertical line at x = 0)
    let ax = 1.0 - smoothstep(0.0, threshold, px.x);
    // ay: coverage of the X-axis (horizontal line at y = 0)
    let ay = 1.0 - smoothstep(0.0, threshold, px.y);

    if ax > ay && ax > 0.01 {
        // Y-axis: muted green — conventional "up" axis in 2-D science plots
        return vec4<f32>(0.25, 0.52, 0.32, ax * 0.85);
    }
    if ay > 0.01 {
        // X-axis: muted red — conventional horizontal axis
        return vec4<f32>(0.55, 0.25, 0.25, ay * 0.85);
    }
    return vec4<f32>(0.0, 0.0, 0.0, 0.0);
}

// ── Fragment ──────────────────────────────────────────────────────────────── //

@fragment
fn fs_grid(@builtin(position) frag_coord : vec4<f32>) -> @location(0) vec4<f32> {

    // World-space coordinate for this fragment.
    let world = (frag_coord.xy - u.center) / u.scale;

    // ── LOD: choose a "nice" fine-grid spacing so lines stay ~60 px apart ── //
    //
    // Using log(x)/log(10) instead of log2(x)/log2(10) avoids the
    // floating-point boundary errors that the mixed-base form can produce
    // near exact powers of ten.
    let raw        = 60.0 / u.scale;
    let log10_raw  = log(abs(raw) + 1e-30) / log(10.0);
    let e          = floor(log10_raw);
    let base       = pow(10.0, e);
    let frac       = raw / base;

    var fine : f32;
    if      frac < 1.5 { fine = base; }
    else if frac < 3.5 { fine = base * 2.0; }
    else if frac < 7.5 { fine = base * 5.0; }
    else               { fine = base * 10.0; }

    let med    = fine * 5.0;   // one medium line every 5 fine lines
    let coarse = fine * 25.0;  // broad orientation lines

    // ── Fine-grid fade ───────────────────────────────────────────────────── //
    // Fine lines fade in only once they are at least ~20 px apart, so they
    // never clutter the view when zoomed far out.
    let fine_px   = fine * u.scale;
    let fine_fade = smoothstep(14.0, 40.0, fine_px);

    // ── Per-level coverage ───────────────────────────────────────────────── //
    let gl_fine   = grid_line(world, fine)   * fine_fade;
    let gl_med    = grid_line(world, med);
    let gl_coarse = grid_line(world, coarse);

    // ── Compositing: brightest level wins ────────────────────────────────── //
    // Each level has its own neutral-blue-gray tint. More prominent levels
    // are progressively brighter. The axis colours (red / green) override
    // the neutral tints when the fragment is on a principal axis.
    var best_a   = 0.0;
    var best_rgb = vec3<f32>(0.15, 0.15, 0.24);

    let a_fine   = gl_fine   * 0.18;
    let a_med    = gl_med    * 0.42;
    let a_coarse = gl_coarse * 0.55;

    if a_fine > best_a {
        best_a   = a_fine;
        best_rgb = vec3<f32>(0.15, 0.15, 0.24);
    }
    if a_med > best_a {
        best_a   = a_med;
        best_rgb = vec3<f32>(0.22, 0.22, 0.32);
    }
    if a_coarse > best_a {
        best_a   = a_coarse;
        best_rgb = vec3<f32>(0.30, 0.30, 0.42);
    }

    // Principal axes override grid color entirely when present.
    let ac = axis_color(world);
    if ac.a > best_a {
        best_a   = ac.a;
        best_rgb = ac.rgb;
    }

    if best_a < 0.004 { discard; }

    return vec4<f32>(best_rgb, best_a);
}

"#;
