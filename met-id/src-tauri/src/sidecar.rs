use serde::Serialize;
use std::fmt;
use tauri::api::process::{Command, CommandEvent};

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
pub fn sidecar_function(sidecar_name: String, sidecar_arguments: Vec<String>) -> Result<String, CommandError> {
    let (mut rx, _child) = Command::new_sidecar(sidecar_name)
        .expect("failed to create `my-sidecar` binary command")
        .args(sidecar_arguments)
        .spawn()
        .expect("Failed to spawn sidecar");
    
    let (tx, rx_data) = std::sync::mpsc::channel();

    tauri::async_runtime::spawn(async move {
        // read events such as stdout
        while let Some(event) = rx.recv().await {
            if let CommandEvent::Stdout(line) = event {
                tx.send(line).unwrap();
            }
        }
    });
    // Here, we collect all the stdout lines. We could also just return the first one, depending on the requirement.
    let mut output = String::new();
    while let Ok(line) = rx_data.recv() {
        output.push_str(&line);
    }


    Ok(output)
}
