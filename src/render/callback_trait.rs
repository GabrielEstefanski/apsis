use eframe::egui_wgpu::{CallbackResources, CallbackTrait};
use std::sync::{Arc, Mutex};

use crate::render::WgpuBackend;

pub struct CallbackFn {
    pub backend: Arc<Mutex<WgpuBackend>>,
    pub device: Arc<wgpu::Device>,
    pub queue: Arc<wgpu::Queue>,
    pub format: wgpu::TextureFormat,
    /// Canvas dimensions in logical pixels: [width, height].
    pub screen: [f32; 2],
    /// Canvas origin in logical pixels: [rect.min.x, rect.min.y].
    pub viewport_min: [f32; 2],
}

impl CallbackTrait for CallbackFn {
    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        pass: &mut wgpu::RenderPass<'static>,
        _resources: &CallbackResources,
    ) {
        let mut backend = self.backend.lock().unwrap();
        let center = backend.center;
        let scale = backend.scale;
        backend.render_frame(
            &self.device,
            &self.queue,
            pass,
            self.screen,
            self.viewport_min,
            self.format,
            center,
            scale,
        );
    }
}
