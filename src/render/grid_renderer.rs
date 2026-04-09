use bytemuck::{Pod, Zeroable};
use std::mem::size_of;

// ── Uniform ───────────────────────────────────────────────────────────────────

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct GridUniform {
    /// Tamanho da tela em pixels.
    screen_size: [f32; 2],
    /// Posição da origem do mundo na tela (pixels desde top-left).
    center: [f32; 2],
    /// Pixels por unidade de mundo.
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

        Self { pipeline, uniform_buf, bind_group }
    }

    pub fn upload(
        &self,
        queue: &wgpu::Queue,
        center: [f32; 2],
        scale: f32,
        screen: [f32; 2],
    ) {
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

// ── Vertex: full-screen quad, nenhum vertex buffer ────────────────────────── //

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

// ── Helpers ──────────────────────────────────────────────────────────────── //

/// Cobertura antialiased de uma linha de grade na coordenada `v` com
/// espaçamento `spacing`. Retorna 0..1 (1 = centro da linha).
fn grid_aa(v : f32, spacing : f32) -> f32 {
    let coord  = v / spacing;
    let deriv  = fwidth(coord);
    let dist   = abs(fract(coord + 0.5) - 0.5);
    return 1.0 - smoothstep(0.0, deriv, dist);
}

/// Máximo entre as linhas X e Y para uma grade 2D.
fn grid_line(world : vec2<f32>, spacing : f32) -> f32 {
    return max(grid_aa(world.x, spacing), grid_aa(world.y, spacing));
}

/// Linha de eixo (x=0 ou y=0), mais grossa.
fn axis_line(world : vec2<f32>) -> f32 {
    let px_dist = abs(world) * u.scale;      // distância em pixels até o eixo
    let threshold = 1.2;                     // meia-largura em pixels
    let ax = 1.0 - smoothstep(0.0, threshold, px_dist.x);
    let ay = 1.0 - smoothstep(0.0, threshold, px_dist.y);
    return max(ax, ay);
}

// ── Fragment ─────────────────────────────────────────────────────────────── //

@fragment
fn fs_grid(@builtin(position) frag_coord : vec4<f32>) -> @location(0) vec4<f32> {

    // Coordenada de mundo correspondente a este pixel
    let world = (frag_coord.xy - u.center) / u.scale;

    // ── LOD: espaçamento "fino" baseado no zoom ───────────────────────────── //
    // Queremos ~60px entre linhas finas → spacing_world = 60 / scale
    let raw = 60.0 / u.scale;

    // Arredonda para número "bonito": 1, 2, 5, 10, 20, 50, ...
    let e     = floor(log2(abs(raw + 1e-30)) / log2(10.0));
    let base  = pow(10.0, e);
    let frac  = raw / base;

    var fine : f32;
    if frac < 1.5       { fine = base; }
    else if frac < 3.5  { fine = base * 2.0; }
    else if frac < 7.5  { fine = base * 5.0; }
    else                { fine = base * 10.0; }

    let med    = fine * 5.0;   // subdivisão maior (1 a cada 5 linhas finas)
    let coarse = fine * 25.0;  // orientação ampla

    // ── Fade do nível fino ──────────────────────────────────────────────────
    // As linhas finas ficam visíveis somente quando >= 20px de distância
    let fine_px   = fine * u.scale;
    let fine_fade = smoothstep(14.0, 40.0, fine_px);

    // ── Cobertura de cada nível ─────────────────────────────────────────────
    let gl_fine   = grid_line(world, fine)   * fine_fade;
    let gl_med    = grid_line(world, med);
    let gl_coarse = grid_line(world, coarse);
    let gl_axis   = axis_line(world);

    // ── Composição: nível mais brilhante vence (evita over-bright) ──────────
    // Eixos: branco/azulado frio  0.44, 0.44, 0.60
    // Coarse: cinza médio         0.30, 0.30, 0.42
    // Med: cinza suave            0.22, 0.22, 0.32
    // Fine: quase invisível       0.15, 0.15, 0.24

    var best_a   = 0.0;
    var best_rgb = vec3<f32>(0.15, 0.15, 0.24);

    let a_fine   = gl_fine   * 0.18;
    let a_med    = gl_med    * 0.42;
    let a_coarse = gl_coarse * 0.55;
    let a_axis   = gl_axis   * 0.78;

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
    if a_axis > best_a {
        best_a   = a_axis;
        best_rgb = vec3<f32>(0.44, 0.44, 0.62);
    }

    if best_a < 0.004 { discard; }

    return vec4<f32>(best_rgb, best_a);
}

"#;
