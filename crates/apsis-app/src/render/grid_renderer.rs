//! Ecliptic-plane grid drawn in world space, with Floating Origin.
//!
//! Full-screen quad whose fragment shader unprojects each pixel
//! through the inverse render-frame view-projection, intersects the
//! camera ray with the absolute-world `z = 0` plane (which sits at
//! `z = -render_origin.z` in the render frame), and shades grid lines
//! on the recovered absolute `(x, y)` coordinates. Same LOD logic as
//! the legacy 2D grid (line spacing chosen so visible lines stay
//! ~60 px apart) but anchored to the world plane the bodies live on
//! rather than to screen-space coordinates.

use bytemuck::{Pod, Zeroable};
use std::mem::size_of;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct GridUniform {
    /// Inverse of the render-frame view-projection. Camera sits at
    /// the origin in this frame.
    inv_view_proj_relative: [[f32; 4]; 4],
    /// Floating-Origin anchor in absolute world units. The ray-plane
    /// intersection runs in the render frame and we add this back to
    /// recover absolute `(x, y)` for the grid LOD math.
    render_origin: [f32; 3],
    _pad0: f32,
    screen_size: [f32; 2],
    _pad1: [f32; 2],
}

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

        Self { pipeline, uniform_buf, bind_group }
    }

    pub fn upload(
        &self,
        queue: &wgpu::Queue,
        inv_view_proj_relative: [[f32; 4]; 4],
        render_origin: [f32; 3],
        screen: [f32; 2],
    ) {
        let u = GridUniform {
            inv_view_proj_relative,
            render_origin,
            _pad0: 0.0,
            screen_size: screen,
            _pad1: [0.0; 2],
        };
        queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&u));
    }

    pub fn draw(&self, pass: &mut wgpu::RenderPass<'_>) {
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.draw(0..6, 0..1);
    }
}

// ── WGSL ─────────────────────────────────────────────────────────────────────

const GRID_SHADER: &str = r#"

struct GridUniform {
    inv_view_proj_relative : mat4x4<f32>,
    render_origin          : vec3<f32>,
    _pad0                  : f32,
    screen_size            : vec2<f32>,
    _pad1                  : vec2<f32>,
};

@group(0) @binding(0)
var<uniform> u : GridUniform;

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

// Antialiased grid line at `v` with spacing `s`. Uses screen-space
// derivative so line width is ~1 px regardless of zoom.
fn grid_aa(v : f32, s : f32) -> f32 {
    let c    = v / s;
    let d    = fwidth(c);
    let dist = abs(fract(c + 0.5) - 0.5);
    return 1.0 - smoothstep(0.0, d, dist);
}

fn grid_line(world : vec2<f32>, s : f32) -> f32 {
    return max(grid_aa(world.x, s), grid_aa(world.y, s));
}

@fragment
fn fs_grid(@builtin(position) frag_coord : vec4<f32>) -> @location(0) vec4<f32> {
    // ── Reconstruct camera ray from this fragment ────────────────────── //
    // Camera sits at the render-frame origin under Floating Origin, so
    // the far-plane unproject *is* the ray direction.
    let ndc = vec2<f32>(
         (frag_coord.x / u.screen_size.x) * 2.0 - 1.0,
        -(frag_coord.y / u.screen_size.y) * 2.0 + 1.0,
    );
    let far_clip     = vec4<f32>(ndc, 0.0, 1.0);
    let far_relative = u.inv_view_proj_relative * far_clip;
    let far_pos      = far_relative.xyz / far_relative.w;
    let ray_dir      = normalize(far_pos);

    // ── Intersect with the absolute z = 0 plane ──────────────────────── //
    // The plane sits at `z_relative = -render_origin.z` in the render
    // frame; reject grazing or upside-down rays so the grid doesn't
    // spit out numerical noise at the horizon.
    if abs(ray_dir.z) < 1e-4 { discard; }
    let target_z = -u.render_origin.z;
    let t = target_z / ray_dir.z;
    if t <= 0.0 { discard; }

    // Hit point in the render frame; recover absolute `(x, y)` for the
    // grid LOD math by adding back the origin. The z component is
    // exactly zero by construction.
    let hit_relative = ray_dir * t;
    let world_xy     = hit_relative.xy + u.render_origin.xy;

    // ── LOD: nice line spacing keyed off screen-space derivative ────── //
    // fwidth on world coordinates gives the "size of one fragment" in
    // world units at this point on the plane — same idea as the 2D
    // grid's `60 / scale`, but generalised to perspective.
    let pixel_size = max(fwidth(world_xy.x), fwidth(world_xy.y));
    let raw        = pixel_size * 60.0;
    let log10_raw  = log(abs(raw) + 1e-30) / log(10.0);
    let e          = floor(log10_raw);
    let base       = pow(10.0, e);
    let frac       = raw / base;

    var fine : f32;
    if      frac < 1.5 { fine = base; }
    else if frac < 3.5 { fine = base * 2.0; }
    else if frac < 7.5 { fine = base * 5.0; }
    else               { fine = base * 10.0; }

    let med    = fine * 5.0;
    let coarse = fine * 25.0;

    // Per-LOD pixel-spacing fade. A LOD fades out before its lines
    // are dense enough to alias to a uniform tint — the same idea
    // applied to fine, med and coarse so coarse lines don't keep
    // painting alpha after they crowd to ~1 px apart.
    let inv_px      = 1.0 / max(pixel_size, 1e-30);
    let fine_px     = fine   * inv_px;
    let med_px      = med    * inv_px;
    let coarse_px   = coarse * inv_px;
    let fine_fade   = smoothstep(14.0, 40.0, fine_px);
    let med_fade    = smoothstep(14.0, 40.0, med_px);
    let coarse_fade = smoothstep(14.0, 40.0, coarse_px);

    let gl_fine   = grid_line(world_xy, fine)   * fine_fade;
    let gl_med    = grid_line(world_xy, med)    * med_fade;
    let gl_coarse = grid_line(world_xy, coarse) * coarse_fade;

    let neutral = vec3<f32>(0.32, 0.34, 0.38);

    let a_fine   = gl_fine   * 0.10;
    let a_med    = gl_med    * 0.22;
    let a_coarse = gl_coarse * 0.32;

    let best_a = max(a_coarse, max(a_med, a_fine));

    // Tight distance fade keyed off camera height + planar offset so
    // the grid stays a local plate of reference instead of bleeding to
    // the horizon. 10× / 30× of the camera's distance from the world
    // origin covers the working volume without painting the whole
    // canvas. `render_origin` is the camera's absolute position, so
    // these expressions match the legacy meaning.
    let r          = length(world_xy);
    let cam_r      = length(u.render_origin.xy);
    let cam_h      = max(cam_r + abs(u.render_origin.z), 1e-3);
    let fade_inner = 10.0 * cam_h;
    let fade_outer = 30.0 * cam_h;
    let dist_fade  = 1.0 - smoothstep(fade_inner, fade_outer, r);

    let a = best_a * dist_fade;
    if a < 0.02 { discard; }

    return vec4<f32>(neutral, a);
}

"#;

#[cfg(test)]
mod shader_tests {
    #[test]
    fn grid_shader_validates() {
        crate::render::validate_wgsl("grid", super::GRID_SHADER);
    }

    /// mat4x4 (64) + vec3 (12) + f32 (4) + vec2 (8) + vec2 (8) = 96 B,
    /// struct alignment 16 (mat4x4 + vec3) → 96 already aligned.
    #[test]
    fn grid_uniform_layout() {
        crate::render::assert_uniform_layout::<super::GridUniform>("GridUniform", 96);
    }
}
