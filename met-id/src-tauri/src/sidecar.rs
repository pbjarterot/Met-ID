use crate::get_app_handle;
use log::warn;
use serde::Serialize;
use std::fmt;
use tauri_plugin_shell::process::CommandEvent;
use tauri_plugin_shell::ShellExt;

// Custom error type
#[derive(Debug, Serialize)]
pub struct CommandError {
    message: String,
}

impl fmt::Display for CommandError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

#[tauri::command]
pub fn sidecar_function(
    sidecar_name: String,
    sidecar_arguments: Vec<String>,
) -> Result<String, CommandError> {
    println!("Starting processing1...");
    let (mut rx, _child) = get_app_handle()
        .unwrap()
        .shell()
        .sidecar(sidecar_name)
        .expect("failed to create `my-sidecar` binary command")
        .args(sidecar_arguments)
        .spawn()
        .expect("Failed to spawn sidecar");

    let (tx, rx_data) = std::sync::mpsc::channel();
    println!("{:?}, {:?}", tx, rx_data);
    tauri::async_runtime::spawn(async move {
        // read events such as stdout
        while let Some(event) = rx.recv().await {
            if let CommandEvent::Stdout(line) = event {
                println!("{:?}", &line);
                tx.send(line).unwrap();
            }
        }
    });
    // Here, we collect all the stdout lines. We could also just return the first one, depending on the requirement.
    let mut output: String = String::new();
    while let Ok(line) = rx_data.recv() {
        println!("{:?}", &line);
        output.push_str(std::str::from_utf8(&line).unwrap());
    }

    Ok(output)
}

pub fn sidecar_function3(
    app: &tauri::AppHandle,
    progress_sender: std::sync::mpsc::Sender<f32>,
    sidecar_arguments: Vec<String>,
) -> std::io::Result<()> {
    println!("Spawning sidecar...");

    let sidecar_command = app
        .shell()
        .sidecar("metabolite")
        .unwrap()
        .args(sidecar_arguments);

    let (mut rx, _child) = sidecar_command.spawn().unwrap();

    tauri::async_runtime::spawn(async move {
        while let Some(event) = rx.recv().await {
            match event {
                CommandEvent::Stdout(line_bytes) => {
                    let line = String::from_utf8_lossy(&line_bytes);
                    println!("Output: {}", line);

                    if let Err(e) = progress_sender.send(100.0) {
                        eprintln!("Failed to send progress: {}", e);
                        warn!("Failed to send progress: {}", e);
                    }
                }
                CommandEvent::Stderr(line_bytes) => {
                    let line = String::from_utf8_lossy(&line_bytes);
                    eprintln!("Error: {}", line);
                    warn!("Error: {}", line);
                }
                _ => {}
            }
        }

        println!("Sidecar finished processing.");
    });

    Ok(())
}
