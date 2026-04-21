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

/// Offscreen colour target for linear-space scene rendering.
pub struct HdrTarget {
    texture: wgpu::Texture,
    view: wgpu::TextureView,
    size: [u32; 2],
    /// Bumps on every reallocation so consumers (e.g. the tonemap bind group)
    /// can detect when their cached `TextureView` is stale.
    generation: u64,
}

impl HdrTarget {
    pub fn new(device: &wgpu::Device, size: [u32; 2]) -> Self {
        let (texture, view) = create_texture(device, size);
        Self { texture, view, size, generation: 1 }
    }

    /// Reallocates the underlying texture when `size` differs from the
    /// currently allocated extent. No-op otherwise. Callers should invoke
    /// this once per frame before recording the scene pass.
    pub fn ensure_size(&mut self, device: &wgpu::Device, size: [u32; 2]) {
        let clamped = [size[0].max(1), size[1].max(1)];
        if clamped == self.size {
            return;
        }
        let (texture, view) = create_texture(device, clamped);
        self.texture = texture;
        self.view = view;
        self.size = clamped;
        self.generation = self.generation.wrapping_add(1);
    }

    #[inline]
    pub fn view(&self) -> &wgpu::TextureView {
        &self.view
    }

    #[inline]
    pub fn size(&self) -> [u32; 2] {
        self.size
    }

    /// Increments on every underlying reallocation. Consumers cache this
    /// value and compare it on each frame to decide when to refresh their
    /// own bindings (e.g. the tonemap pipeline's sampled-texture bind group).
    #[inline]
    pub fn generation(&self) -> u64 {
        self.generation
    }
}

fn create_texture(device: &wgpu::Device, size: [u32; 2]) -> (wgpu::Texture, wgpu::TextureView) {
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
