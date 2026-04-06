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
            let x = physical_px.max(0.0);

            let log = (x + 1.0).ln();
            let gamma = 0.6;
            let scaled = log.powf(gamma) * 10.0;

            scaled.max(params.min_px).max(2.5)
        }
    }
}
