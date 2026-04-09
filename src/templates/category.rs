use eframe::egui;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TemplateCategory {
    Bodies,
    Systems,
    ThreeBodyProblems,
}

impl TemplateCategory {
    pub fn label(self) -> &'static str {
        match self {
            Self::Bodies => "BODIES",
            Self::Systems => "SYSTEMS",
            Self::ThreeBodyProblems => "3-BODY PROBLEMS",
        }
    }

    pub fn grid_id(self) -> egui::Id {
        egui::Id::new(format!("cat_grid_{:?}", self))
    }
}
