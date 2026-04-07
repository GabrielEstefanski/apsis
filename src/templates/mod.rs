pub mod builders;
pub mod catalog;
pub mod category;
pub mod core;
pub mod instantiate;
pub mod presets;

pub use catalog::{TEMPLATES, TemplateEntry};
pub use category::TemplateCategory;
pub use core::{Template, TemplateBody};
pub use instantiate::{instantiate, instantiate_at};
