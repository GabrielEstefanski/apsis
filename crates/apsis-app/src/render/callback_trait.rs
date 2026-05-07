use eframe::egui_wgpu::{CallbackResources, CallbackTrait, ScreenDescriptor};
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
    /// World → clip transform driving the body pass. Column-major.
    pub view_proj: [[f32; 4]; 4],
    /// Eye position in world space — fragment shader needs it for
    /// ray-sphere intersection.
    pub camera_pos: [f32; 3],
}

impl CallbackTrait for CallbackFn {
    /// Records the scene into the HDR offscreen target **before** egui
    /// begins its main render pass. The supplied `egui_encoder` is submitted
    /// ahead of the swapchain pass, so anything we record here lands in the
    /// HDR texture and is available for the tonemap composite in `paint`.
    fn prepare(
        &self,
        _device: &wgpu::Device,
        queue: &wgpu::Queue,
        screen_descriptor: &ScreenDescriptor,
        egui_encoder: &mut wgpu::CommandEncoder,
        _callback_resources: &mut CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        let ppp = screen_descriptor.pixels_per_point;
        let physical_size = [
            ((self.screen[0] * ppp).ceil() as u32).max(1),
            ((self.screen[1] * ppp).ceil() as u32).max(1),
        ];

        let mut backend = self.backend.lock().unwrap();
        backend.prepare_scene(
            &self.device,
            queue,
            egui_encoder,
            self.screen,
            self.viewport_min,
            physical_size,
            ppp,
            self.format,
            self.view_proj,
            self.camera_pos,
        );
        Vec::new()
    }

    /// Tonemap composite: samples the HDR target recorded in `prepare` and
    /// writes the result into the canvas region of the swapchain attachment.
    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        pass: &mut wgpu::RenderPass<'static>,
        _resources: &CallbackResources,
    ) {
        let backend = self.backend.lock().unwrap();
        backend.composite(pass);
    }
}
