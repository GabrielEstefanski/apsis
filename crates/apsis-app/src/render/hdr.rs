//! HDR offscreen targets — intermediate buffers for linear-space rendering.
//!
//! Two `Rgba16Float` colour textures share a single `Depth32Float`
//! attachment. The **reflective** plane (`view_r`) carries grid, trails,
//! lines, circles, point sprites, and the lit hemisphere of every body;
//! the auto-exposure metering reads it. The **luminous** plane (`view_l`)
//! carries the emissive contribution of self-luminous bodies and feeds
//! the bloom pass.
//!
//! The composite pass samples both, multiplies the reflective half by the
//! exposure scalar, sums into one HDR signal, and runs a single ACES
//! curve.
//!
//! # Format
//!
//! `Rgba16Float` keeps ~11 bits of mantissa across a wide dynamic range —
//! enough headroom for additive light accumulation without clipping.
//!
//! # Ownership
//!
//! The targets are owned by [`crate::render::WgpuBackend`] and grow/shrink
//! together via [`HdrTarget::ensure_size`]. A single generation counter
//! covers both colour views and the depth view.

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

/// Offscreen render targets for linear-space scene rendering.
///
/// Holds the reflective and luminous colour planes plus the shared depth
/// attachment. Reallocation is atomic across all three textures and bumps
/// `generation` so downstream consumers know to rebuild their bind groups.
pub struct HdrTarget {
    color_r_texture: wgpu::Texture,
    color_r_view: wgpu::TextureView,
    color_l_texture: wgpu::Texture,
    color_l_view: wgpu::TextureView,
    depth_texture: wgpu::Texture,
    depth_view: wgpu::TextureView,
    size: [u32; 2],
    generation: u64,
}

impl HdrTarget {
    pub fn new(device: &wgpu::Device, size: [u32; 2]) -> Self {
        let (color_r_texture, color_r_view) = create_color(device, size, "hdr::scene_color_r");
        let (color_l_texture, color_l_view) = create_color(device, size, "hdr::scene_color_l");
        let (depth_texture, depth_view) = create_depth(device, size);
        Self {
            color_r_texture,
            color_r_view,
            color_l_texture,
            color_l_view,
            depth_texture,
            depth_view,
            size,
            generation: 1,
        }
    }

    /// Reallocates the underlying textures when `size` differs from the
    /// currently allocated extent. No-op otherwise.
    pub fn ensure_size(&mut self, device: &wgpu::Device, size: [u32; 2]) {
        let clamped = [size[0].max(1), size[1].max(1)];
        if clamped == self.size {
            return;
        }
        let (color_r_texture, color_r_view) = create_color(device, clamped, "hdr::scene_color_r");
        let (color_l_texture, color_l_view) = create_color(device, clamped, "hdr::scene_color_l");
        let (depth_texture, depth_view) = create_depth(device, clamped);
        self.color_r_texture = color_r_texture;
        self.color_r_view = color_r_view;
        self.color_l_texture = color_l_texture;
        self.color_l_view = color_l_view;
        self.depth_texture = depth_texture;
        self.depth_view = depth_view;
        self.size = clamped;
        self.generation = self.generation.wrapping_add(1);
    }

    /// Reflective colour plane. Sampled by the metering reducer and the
    /// composite pass.
    #[inline]
    pub fn view_r(&self) -> &wgpu::TextureView {
        &self.color_r_view
    }

    /// Luminous colour plane. Sampled by the bloom and composite passes.
    #[inline]
    pub fn view_l(&self) -> &wgpu::TextureView {
        &self.color_l_view
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

fn create_color(
    device: &wgpu::Device,
    size: [u32; 2],
    label: &str,
) -> (wgpu::Texture, wgpu::TextureView) {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some(label),
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
