use std::io::{self, BufRead};
use std::sync::mpsc::{self, Receiver};
use std::thread;

pub fn spawn_input_thread() -> Receiver<String> {
    let (tx, rx) = mpsc::channel();

    thread::spawn(move || {
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            if let Ok(cmd) = line {
                tx.send(cmd).ok();
            }
        }
    });

    rx
}