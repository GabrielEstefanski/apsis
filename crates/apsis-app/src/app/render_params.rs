use crate::app::ui::SemanticScaleMode;

#[derive(Clone, Copy)]
pub struct RenderParams {
    pub world_scale: f32,
    pub mode: SemanticScaleMode,
    pub min_px: f32,
}

pub fn compute_render_radius(physical_radius: f64, params: RenderParams) -> f32 {
    let physical_px = physical_radius as f32 * params.world_scale;

    match params.mode {
        SemanticScaleMode::Physical => physical_px,

        SemanticScaleMode::Comparative => physical_px.max(params.min_px),

        SemanticScaleMode::Illustrative => {
            let k = 0.15;
            let scaled = (1.0 - (-k * physical_px).exp()) * 20.0;

            scaled.max(params.min_px).max(2.5)
        },
    }
}
