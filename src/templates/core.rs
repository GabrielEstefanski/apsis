use crate::domain::materials::Material;

#[derive(Clone, Copy)]
pub struct TemplateBody {
    pub mass: f64,
    pub radius: f64,
    pub position: Option<[f64; 2]>,
    pub velocity: [f64; 2],
    pub material: Material,
}

pub struct Template {
    pub name: &'static str,
    pub bodies: Vec<TemplateBody>,

    pub scale: f64,
}
