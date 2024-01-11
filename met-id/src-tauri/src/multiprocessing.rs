use std::thread;
use std::sync::mpsc;
use std::collections::HashMap;
use crate::add_to_db::add_to_db_functions::add_to_db_rust;
use tauri::Window;

#[derive(serde::Deserialize)]
#[allow(non_snake_case)]
pub struct MyCommandArgs {
    name: String, 
    smilesSmartsMz: String, 
    metType: String, 
    endoExoOrOther: HashMap<String, bool>, 
    inTissue: HashMap<String, bool>, 
    adducts: HashMap<String, String>
}

#[tauri::command]
pub fn my_command(window: Window, args: MyCommandArgs) {
    let (tx, rx) = mpsc::channel();
    println!("my_command");
    thread::spawn(move || {
        add_to_db_rust(args.name, args.smilesSmartsMz, args.metType, args.endoExoOrOther, args.inTissue, args.adducts, &tx);
        // Signal completion using -1.0 as the sentinel value
        tx.send(-1.0).unwrap();
    });

    // Spawn a new thread for listening to progress updates
    thread::spawn(move || {
        loop {
            match rx.recv() {
                Ok(progress) => {
                    if (progress + 1.0).abs() < std::f32::EPSILON { // Check if progress is close to -1.0
                        break;  // Exit the loop once the task is done
                    }
                    // Update the frontend using the `window` object.
                    window.emit("update-progress", Some(&progress.to_string())).unwrap();
                },
                Err(e) => {
                    eprintln!("Error receiving progress: {}", e);
                    break;
                }
            }
        }
    });
}
