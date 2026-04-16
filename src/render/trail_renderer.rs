use bytemuck::{Pod, Zeroable};
use std::mem::size_of;

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
    pub decay_k: f32,
    pub tail_desaturate: f32,
    pub _pad: [f32; 2],
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
    uploaded_sample_count: u64,
    uploaded_len: u32,
    visible_colors: Vec<[f32; 4]>,
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
                wgpu::BindGroupEntry { binding: 0, resource: pos_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: color_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: state_buf.as_entire_binding() },
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
            uploaded_sample_count: 0,
            uploaded_len: 0,
            visible_colors: Vec::new(),
        }
    }

    pub fn upload(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        trail: &TrailBuffer,
        visibility: Option<&[bool]>,
        center: [f32; 2],
        scale: f32,
        trail_width: f32,
        decay_k: f32,         // novo
        tail_desaturate: f32, // novo
    ) {
        let n_bodies = trail.n_bodies();
        let cap = trail.capacity();
        let sample_count = trail.sample_count();

        if self.n_bodies != n_bodies || self.cap != cap {
            self.n_bodies = n_bodies;
            self.cap = cap;
            self.uploaded_sample_count = 0;
            self.uploaded_len = 0;

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
                    wgpu::BindGroupEntry { binding: 0, resource: self.pos_buf.as_entire_binding() },
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

        let delta = sample_count.saturating_sub(self.uploaded_sample_count);
        let must_full_upload = self.uploaded_sample_count == 0
            || delta >= cap as u64
            || trail.len() < self.uploaded_len;

        if must_full_upload {
            queue.write_buffer(&self.pos_buf, 0, bytemuck::cast_slice(trail.positions()));
        } else if delta > 0 {
            let col_bytes = n_bodies as u64 * size_of::<[f32; 2]>() as u64;
            let first_col = (trail.head() + cap - delta as u32) % cap;

            for step in 0..delta as u32 {
                let col = (first_col + step) % cap;
                let offset = col as u64 * col_bytes;
                queue.write_buffer(
                    &self.pos_buf,
                    offset,
                    bytemuck::cast_slice(trail.column_slice(col)),
                );
            }
        }

        if let Some(visibility) = visibility {
            self.visible_colors.clear();
            self.visible_colors.extend(trail.colors().iter().enumerate().map(|(i, color)| {
                let mut rgba = *color;
                if !visibility.get(i).copied().unwrap_or(true) {
                    rgba[3] = 0.0;
                }
                rgba
            }));
            queue.write_buffer(&self.color_buf, 0, bytemuck::cast_slice(&self.visible_colors));
        } else {
            queue.write_buffer(&self.color_buf, 0, bytemuck::cast_slice(trail.colors()));
        }

        let state = TrailState {
            head: trail.head(),
            trail_len: trail.len(),
            n_bodies,
            cap,
            center,
            scale,
            trail_width,
            decay_k,
            tail_desaturate,
            _pad: [0.0; 2],
        };

        queue.write_buffer(&self.state_buf, 0, bytemuck::bytes_of(&state));
        self.uploaded_sample_count = sample_count;
        self.uploaded_len = trail.len();
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
    viewport_min: vec2<f32>,
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
    decay_k: f32,
    tail_desaturate: f32,
    _pad: vec2<f32>,
};

@group(1) @binding(0) var<storage, read> positions: array<vec2<f32>>;
@group(1) @binding(1) var<storage, read> colors: array<vec4<f32>>;
@group(1) @binding(2) var<uniform> state: TrailState;

// ── Utilitários ────────────────────────────────────────────────────────────────

fn to_ndc(p: vec2<f32>) -> vec4<f32> {
    let local = p - screen.viewport_min;
    let x =  (local.x / screen.size.x) * 2.0 - 1.0;
    let y = -(local.y / screen.size.y) * 2.0 + 1.0;
    return vec4<f32>(x, y, 0.0, 1.0);
}

// Converte RGB → luminância (BT.709)
fn luminance(rgb: vec3<f32>) -> f32 {
    return dot(rgb, vec3<f32>(0.2126, 0.7152, 0.0722));
}

// Desatura parcialmente uma cor: t=0 → original, t=1 → tons de cinza
fn desaturate(rgb: vec3<f32>, t: f32) -> vec3<f32> {
    let lum = luminance(rgb);
    return mix(rgb, vec3<f32>(lum), t);
}

// Shift de temperatura: age=0 (cauda) → frio/escuro, age=1 (ponta) → quente/brilhante
// Simula emissão física como Universe Sandbox
fn temperature_tint(rgb: vec3<f32>, age: f32) -> vec3<f32> {
    // Ponta (age≈1): ligeiramente mais quente → boost no vermelho/verde, reduz azul
    // Cauda (age≈0): esfria e desatura
    let warm  = vec3<f32>(1.12, 1.05, 0.82);  // tint quente
    let cold  = vec3<f32>(0.72, 0.78, 0.92);  // tint frio/azulado
    let tint  = mix(cold, warm, age);

    // Boost de saturação/brilho na ponta
    let boosted = rgb * tint;

    // Desaturação progressiva na cauda
    let desat_amount = (1.0 - age) * state.tail_desaturate;
    return desaturate(boosted, desat_amount);
}

struct VSOut {
    @builtin(position) pos: vec4<f32>,
    @location(0)       color: vec4<f32>,
};

@vertex
fn vs_trail(
    @builtin(vertex_index)   vi:   u32,
    @builtin(instance_index) body: u32,
) -> VSOut {
    let seg = vi / 6u;
    let tri = vi % 6u;

    let i0 = (state.head + seg)      % state.cap;
    let i1 = (state.head + seg + 1u) % state.cap;

    let idx0 = i0 * state.n_bodies + body;
    let idx1 = i1 * state.n_bodies + body;

    let p0 = positions[idx0];
    let p1 = positions[idx1];

    // Descarta segmentos inválidos (NaN ≠ NaN)
    if any(p0 != p0) || any(p1 != p1) {
        var out: VSOut;
        out.pos   = vec4<f32>(0.0);
        out.color = vec4<f32>(0.0);
        return out;
    }

    let screen_p0 = state.center + p0 * state.scale;
    let screen_p1 = state.center + p1 * state.scale;

    let dir     = screen_p1 - screen_p0;
    let len     = max(length(dir), 1e-5);
    let tangent = dir / len;
    let normal  = vec2<f32>(-tangent.y, tangent.x);

    // ── Idade normalizada [0 = cauda/mais antiga → 1 = ponta/mais recente] ──
    let oldest_seg = state.cap - state.trail_len;
    let raw_age    = f32(i32(seg) - i32(oldest_seg))
                   / f32(max(state.trail_len, 2u) - 1u);
    let age = clamp(raw_age, 0.0, 1.0);

    // ── Largura variável: estreita na cauda, larga na ponta ──────────────────
    // Curva pow dá crescimento suave e visualmente convincente
    let width_factor = pow(age, 0.55);                       // 0→0 … 1→1
    let half_width   = state.trail_width * width_factor;

    // ── Quad vertices ────────────────────────────────────────────────────────
    var uv = array<vec2<f32>, 6>(
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0, -1.0),
        vec2<f32>(1.0,  1.0),
        vec2<f32>(0.0,  1.0),
    );
    let corner = uv[tri];
    let pos = screen_p0
            + tangent * (len     * corner.x)
            + normal  * (half_width * corner.y);

    // ── Alpha com decaimento exponencial (estilo simulador físico) ───────────
    //
    // exp(-k * (1 - age)):
    //   age=1 (ponta)  → exp(0)   = 1.0  (opaco)
    //   age=0 (cauda)  → exp(-k)  ≈ 0.0  (transparente)
    //
    // k=5 → cauda some suavemente; k=8 → corte mais agressivo
    // Fade-in curto na cauda para não ter corte duro quando o buffer está cheio
    let decay    = exp(-state.decay_k * (1.0 - age));
    let fade_in  = smoothstep(0.0, 0.04, age);   // apaga os ~4% mais antigos
    let alpha    = colors[body].a * decay * fade_in;

    // ── Cor com temperatura e desaturação na cauda ───────────────────────────
    let base_rgb  = colors[body].rgb;
    let tinted    = temperature_tint(base_rgb, age);

    // Clareia ligeiramente a ponta (simula core quente)
    let core_boost = mix(1.0, 1.25, smoothstep(0.82, 1.0, age));
    let final_rgb  = clamp(tinted * core_boost, vec3<f32>(0.0), vec3<f32>(1.5));

    var out: VSOut;
    out.pos   = to_ndc(pos);
    out.color = vec4<f32>(final_rgb, alpha);
    return out;
}

@fragment
fn fs_trail(in: VSOut) -> @location(0) vec4<f32> {
    if in.color.a <= 0.004 { discard; }
    return in.color;
}
"#;
