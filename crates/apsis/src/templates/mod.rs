pub mod builders;
pub mod catalog;
pub mod category;
pub mod core;
pub mod instantiate;
pub mod keplerian;
pub mod kind;
pub mod presets;

pub use catalog::TEMPLATES;
pub use category::TemplateCategory;
pub use core::{Template, TemplateBody, UnitSystem};
pub use instantiate::instantiate_at;
pub use kind::{TemplateKind, UnknownTemplate};
