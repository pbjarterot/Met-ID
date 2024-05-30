use serde::Serialize;
use std::fmt;
use tauri::api::process::{Command, CommandEvent};
use std::sync::{Arc, Mutex};
use std::thread;
use std::sync::mpsc::{self, Sender, Receiver};
use tauri::async_runtime::block_on;


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
    println!("Starting processing1...");
    let (mut rx, _child) = Command::new_sidecar(sidecar_name)
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
        output.push_str(&line);
    }

    Ok(output)
}

#[tauri::command]
pub fn sidecar_function2(sidecar_name: String, sidecar_arguments: Vec<String>) -> Result<String, CommandError> {
    println!("Starting processing2...");
    let (mut rx, _child) = Command::new_sidecar(sidecar_name)
        .expect("failed to create `my-sidecar` binary command")
        .args(sidecar_arguments)
        .spawn()
        .expect("Failed to spawn sidecar");

    // Channel for sending stdout lines from the async runtime to the main thread
    let (tx, rx_data): (Sender<String>, Receiver<String>) = mpsc::channel();

    // Arc and Mutex to share the sender between threads safely
    let tx = Arc::new(Mutex::new(tx));

    // Spawn a thread to handle the async runtime
    let tx_clone = Arc::clone(&tx);
    thread::spawn(move || {
        block_on(async move {
            while let Some(event) = rx.recv().await {
                if let CommandEvent::Stdout(line) = event {
                    println!("From async thread: {:?}", &line);
                    tx_clone.lock().unwrap().send(line).unwrap();
                }
            }
        });
    });

    // Collect all the stdout lines in the main thread
    let mut output: String = String::new();
    while let Ok(line) = rx_data.recv() {
        println!("From main thread: {:?}", &line); // This is just for debugging, can be removed
        output.push_str(&line);
    }

    Ok(output)
}