use crate::add_to_db::add_to_db_rust;
use log::warn;
use std::collections::HashMap;
use std::sync::mpsc;
use std::thread;
use tauri::{Emitter, Window};

#[derive(serde::Deserialize)]
#[allow(non_snake_case)]
pub struct MyCommandArgs {
    name: String,
    smilesSmartsMz: String,
    metType: String,
    endoExoOrOther: HashMap<String, bool>,
    inTissue: HashMap<String, bool>,
    adducts: Vec<String>,
}

#[tauri::command]
pub fn my_command(window: Window, args: MyCommandArgs) {
    let (tx, rx) = mpsc::channel();
    let tx_clone = tx.clone();
    thread::spawn(move || {
        tx.send(1 as f32).unwrap();
        add_to_db_rust(
            args.name,
            args.smilesSmartsMz,
            args.metType,
            args.endoExoOrOther,
            args.inTissue,
            args.adducts,
            tx,
        );
        // Signal completion using -1.0 as the sentinel value
        tx_clone.send(-1.0).unwrap();
    });

    // Spawn a new thread for listening to progress updates
    thread::spawn(move || {
        loop {
            match rx.recv() {
                Ok(progress) => {
                    if (progress + 1.0).abs() < std::f32::EPSILON {
                        // Check if progress is close to -1.0
                        break; // Exit the loop once the task is done
                    }
                    // Update the frontend using the `window` object.
                    window
                        .emit("update-progress", Some(&progress.to_string()))
                        .unwrap();
                }
                Err(e) => {
                    eprintln!("Error receiving progress: {}", e);
                    warn!("error receiving progress: {}", e);
                    break;
                }
            }
        }
    });
}
