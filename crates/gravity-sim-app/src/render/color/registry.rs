//! Registries for the colour pipeline. Follow the same pattern as
//! [`FieldRegistry`](gravity_sim_core::domain::field::FieldRegistry): own the trait
//! objects, resolve by `id()`.

use super::colormap::Colormap;
use super::normalizer::Normalizer;

pub struct ColormapRegistry {
    entries: Vec<Box<dyn Colormap>>,
}

impl ColormapRegistry {
    pub fn new() -> Self {
        Self { entries: Vec::new() }
    }

    pub fn standard() -> Self {
        use super::cool_warm::CoolWarm;
        use super::grayscale::Grayscale;
        use super::inferno::Inferno;
        use super::plasma::Plasma;
        use super::viridis::Viridis;

        let mut r = Self::new();
        r.register(Box::new(Viridis));
        r.register(Box::new(Inferno));
        r.register(Box::new(Plasma));
        r.register(Box::new(CoolWarm));
        r.register(Box::new(Grayscale));
        r
    }

    pub fn register(&mut self, cm: Box<dyn Colormap>) {
        self.entries.push(cm);
    }

    pub fn get(&self, id: &str) -> Option<&dyn Colormap> {
        self.entries.iter().find(|c| c.id() == id).map(|c| c.as_ref())
    }

    pub fn iter(&self) -> impl Iterator<Item = &dyn Colormap> {
        self.entries.iter().map(|c| c.as_ref())
    }
}

impl Default for ColormapRegistry {
    fn default() -> Self {
        Self::standard()
    }
}

pub struct NormalizerRegistry {
    entries: Vec<Box<dyn Normalizer>>,
}

impl NormalizerRegistry {
    pub fn new() -> Self {
        Self { entries: Vec::new() }
    }

    pub fn standard() -> Self {
        use super::normalizer::{LinearNormalizer, LogNormalizer};

        let mut r = Self::new();
        r.register(Box::new(LinearNormalizer));
        r.register(Box::new(LogNormalizer));
        r
    }

    pub fn register(&mut self, n: Box<dyn Normalizer>) {
        self.entries.push(n);
    }

    pub fn get(&self, id: &str) -> Option<&dyn Normalizer> {
        self.entries.iter().find(|n| n.id() == id).map(|n| n.as_ref())
    }

    pub fn iter(&self) -> impl Iterator<Item = &dyn Normalizer> {
        self.entries.iter().map(|n| n.as_ref())
    }
}

impl Default for NormalizerRegistry {
    fn default() -> Self {
        Self::standard()
    }
}
