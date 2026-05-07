//! HDR offscreen target — intermediate buffer for linear-space rendering.
//!
//! The simulation's scene passes (grid, trails, bodies, orbit overlays) draw
//! into an `Rgba16Float` texture sized to the canvas in physical pixels. A
//! separate tonemap pass (see [`crate::render::tonemap`]) samples this
//! texture, applies exposure + ACES, and outputs to the sRGB swapchain.
//!
//! # Why `Rgba16Float`
//!
//! Linear-space lighting math accumulates values that legitimately exceed 1.0
//! (multi-star scenes, emissive bodies). A `Unorm8` surface clips them and
//! loses all highlight detail. `Rgba16Float` keeps roughly 11 bits of mantissa
//! across a vast dynamic range — enough headroom for additive light
//! accumulation without artefacts.
//!
//! # Ownership
//!
//! The target is owned by [`crate::render::WgpuBackend`] and grows/shrinks
//! with the canvas via [`HdrTarget::ensure_size`]. A generation counter
//! increments on every resize so the tonemap pipeline knows when to rebuild
//! its sampled-texture bind group.

/// Intermediate HDR format for scene rendering.
///
/// Chosen over `Rgba32Float` because 16-bit floats already provide ample
/// headroom for the simulation's brightness range while halving bandwidth.
pub const HDR_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba16Float;

/// Companion depth target for the body pass. `Depth32Float` paired with
/// reverse-Z infinite-far gives essentially uniform precision across the
/// AU-to-light-year range — the standard modern choice for astronomical
/// scales.
pub const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

/// Offscreen colour target for linear-space scene rendering.
///
/// Owns a depth view alongside the colour view so both stay sized together
/// and a single generation counter signals reallocation to consumers.
pub struct HdrTarget {
    color_texture: wgpu::Texture,
    color_view: wgpu::TextureView,
    depth_texture: wgpu::Texture,
    depth_view: wgpu::TextureView,
    size: [u32; 2],
    generation: u64,
}

impl HdrTarget {
    pub fn new(device: &wgpu::Device, size: [u32; 2]) -> Self {
        let (color_texture, color_view) = create_color(device, size);
        let (depth_texture, depth_view) = create_depth(device, size);
        Self { color_texture, color_view, depth_texture, depth_view, size, generation: 1 }
    }

    /// Reallocates the underlying textures when `size` differs from the
    /// currently allocated extent. No-op otherwise.
    pub fn ensure_size(&mut self, device: &wgpu::Device, size: [u32; 2]) {
        let clamped = [size[0].max(1), size[1].max(1)];
        if clamped == self.size {
            return;
        }
        let (color_texture, color_view) = create_color(device, clamped);
        let (depth_texture, depth_view) = create_depth(device, clamped);
        self.color_texture = color_texture;
        self.color_view = color_view;
        self.depth_texture = depth_texture;
        self.depth_view = depth_view;
        self.size = clamped;
        self.generation = self.generation.wrapping_add(1);
    }

    #[inline]
    pub fn view(&self) -> &wgpu::TextureView {
        &self.color_view
    }

    #[inline]
    pub fn depth_view(&self) -> &wgpu::TextureView {
        &self.depth_view
    }

    #[inline]
    pub fn size(&self) -> [u32; 2] {
        self.size
    }

    #[inline]
    pub fn generation(&self) -> u64 {
        self.generation
    }
}

fn create_color(device: &wgpu::Device, size: [u32; 2]) -> (wgpu::Texture, wgpu::TextureView) {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("hdr::scene_color"),
        size: wgpu::Extent3d { width: size[0], height: size[1], depth_or_array_layers: 1 },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: HDR_FORMAT,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });
    let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    (texture, view)
}

fn create_depth(device: &wgpu::Device, size: [u32; 2]) -> (wgpu::Texture, wgpu::TextureView) {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("hdr::scene_depth"),
        size: wgpu::Extent3d { width: size[0], height: size[1], depth_or_array_layers: 1 },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: DEPTH_FORMAT,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    });
    let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    (texture, view)
}
