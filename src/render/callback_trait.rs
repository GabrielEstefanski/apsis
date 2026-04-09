use eframe::egui_wgpu::{CallbackResources, CallbackTrait};
use std::sync::{Arc, Mutex};

use crate::render::WgpuBackend;

pub struct CallbackFn {
    pub backend: Arc<Mutex<WgpuBackend>>,
    pub device: Arc<wgpu::Device>,
    pub queue: Arc<wgpu::Queue>,
    pub format: wgpu::TextureFormat,
    pub screen: [f32; 2],
}

impl CallbackTrait for CallbackFn {
    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        pass: &mut wgpu::RenderPass<'static>,
        _resources: &CallbackResources,
    ) {
        let mut backend = self.backend.lock().unwrap();

        let (trail_ptr, center, scale) = {
            let trail_ptr = backend
                .trail_buffer
                .as_ref()
                .map(|t| t as *const crate::core::trail_buffer::TrailBuffer);

            (trail_ptr, backend.center, backend.scale)
        };

        let trail_buf = unsafe { trail_ptr.map(|p| &*p) };

        backend.render_frame(
            &self.device,
            &self.queue,
            pass,
            self.screen,
            self.format,
            trail_buf,
            center,
            scale,
        );
    }
}
