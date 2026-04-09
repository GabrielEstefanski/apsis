pub mod egui_backend;
pub mod render_backend;
pub mod trail_backend;
pub mod wgpu_backend;

pub use render_backend::RenderBackend;
pub use trail_backend::TrailRenderer;
pub use wgpu_backend::WgpuBackend;
