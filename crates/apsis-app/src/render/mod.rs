pub mod bloom;
pub mod callback_trait;
pub mod color;
pub mod exposure;
pub mod grid_renderer;
pub mod hdr;
pub mod lighting;
pub mod luminance_reducer;
pub mod orbit_overlay;
pub mod orbit_smoother;
pub mod point_renderer;
pub mod render_relative;
pub mod tonemap;
pub mod trail;
pub mod trail_renderer;
pub mod wgpu_backend;

pub use callback_trait::CallbackFn;
pub use exposure::ExposureState;
pub use lighting::{LightSpec, SceneLighting};
pub use trail::{TrailStyle, TrailStylePreset};
pub use trail_renderer::TrailRenderer;
pub use wgpu_backend::WgpuBackend;

/// Parse + validate a WGSL shader source through `naga`.
///
/// `device.create_shader_module` only runs at app startup, so an
/// identifier-scope or type error in a shader string slips past
/// `cargo build` / `cargo clippy` and only surfaces as a runtime
/// `wgpu_core` validation error. Each renderer module is expected to
/// invoke this helper on its shader from a `#[test]` so the same
/// failure modes fail at `cargo test`.
#[cfg(test)]
pub(crate) fn validate_wgsl(name: &str, source: &str) {
    let module = naga::front::wgsl::parse_str(source)
        .unwrap_or_else(|e| panic!("WGSL parse failed for `{name}`:\n{e}"));
    let mut validator = naga::valid::Validator::new(
        naga::valid::ValidationFlags::all(),
        naga::valid::Capabilities::all(),
    );
    validator
        .validate(&module)
        .unwrap_or_else(|e| panic!("WGSL validation failed for `{name}`:\n{e:?}"));
}
