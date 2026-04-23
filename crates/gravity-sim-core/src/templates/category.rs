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
}
