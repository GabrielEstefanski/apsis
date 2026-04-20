pub mod callback_trait;
pub mod grid_renderer;
pub mod trail;
pub mod trail_buffer;
pub mod trail_renderer;
pub mod wgpu_backend;

pub use callback_trait::CallbackFn;
pub use trail::{TrailStyle, TrailStylePreset};
pub use trail_renderer::TrailRenderer;
pub use wgpu_backend::WgpuBackend;
