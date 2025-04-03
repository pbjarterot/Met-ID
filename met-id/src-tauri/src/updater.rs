use tauri_plugin_updater::UpdaterExt;
use tauri::{AppHandle, Emitter};
use std::sync::{Arc, Mutex};
use std::collections::HashMap;
use std::sync::mpsc;
use uuid::Uuid;
use log::info;
// For storing pending responses
pub type CallbackMap = Arc<Mutex<HashMap<String, mpsc::Sender<bool>>>>;

#[tauri::command]
pub fn frontend_bool_response(id: String, value: bool, state: tauri::State<CallbackMap>) {
    if let Some(sender) = state.lock().unwrap().remove(&id) {
        sender.send(value).unwrap();
    }
}

pub fn emit_and_wait(app: &AppHandle, state: CallbackMap) -> bool {
    let (tx, rx) = mpsc::channel();
    let id = Uuid::new_v4().to_string();

    // Store the sender with the unique ID
    state.lock().unwrap().insert(id.clone(), tx);

    // Emit event with the unique ID
    app.emit("ask-frontend-bool", id.clone()).unwrap();

    // Wait for response (blocking)
    let response = rx.recv().unwrap();
    println!("Got response from frontend: {}", response);
    response
}


pub async fn update(app: tauri::AppHandle, callback_map: CallbackMap) -> tauri_plugin_updater::Result<()> {
    if let Some(update) = app.updater()?.check().await? {

        let state = callback_map.clone();
        let update_bool = emit_and_wait(&app, state);
        info!("Result from frontend: {}", update_bool);

        if update_bool == true {
            let mut downloaded = 0;

            // alternatively we could also call update.download() and update.install() separately
            update
            .download_and_install(
                |chunk_length, content_length| {
                downloaded += chunk_length;
                println!("downloaded {downloaded} from {content_length:?}");
                },
                || {
                println!("download finished");
                },
            )
            .await?;

            println!("update installed");
            app.restart();
        }
    }

    Ok(())
}