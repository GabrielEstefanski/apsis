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

/// Asserts that the Rust `#[repr(C)]` size of an instance / uniform
/// struct matches the size the WGSL counterpart will demand at
/// bind-group validation time. Catches the failure mode where a
/// `vec3<f32>` field forces the WGSL struct alignment up to 16 and
/// the Rust side stays packed tighter, producing the runtime error
/// "buffer bound at binding index N is bound with size A where the
/// shader expects B".
///
/// `expected_wgsl_size` is computed by hand from the WGSL alignment
/// rules; the test fails loudly with both numbers if they disagree.
#[cfg(test)]
pub(crate) fn assert_uniform_layout<T>(name: &str, expected_wgsl_size: usize) {
    let rust_size = std::mem::size_of::<T>();
    assert_eq!(
        rust_size, expected_wgsl_size,
        "Rust `#[repr(C)]` size of `{name}` ({rust_size} B) disagrees with the \
         WGSL struct size the shader expects ({expected_wgsl_size} B). Adjust \
         the trailing `_pad` to bring them in line, or rework the WGSL struct."
    );
}
