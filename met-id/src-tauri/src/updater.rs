use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use tauri::{AppHandle, Emitter, State};
use tauri_plugin_updater::UpdaterExt;
//use std::sync::mpsc;
use log::info;
use tokio::sync::oneshot;

// For storing pending responses
pub type CallbackMap = Arc<Mutex<HashMap<String, oneshot::Sender<bool>>>>;

#[tauri::command]
pub fn frontend_bool_response(id: String, value: bool, state: tauri::State<CallbackMap>) {
    println!("[frontend_bool_response] Received callback for ID: {}", id);

    let mut map = state.lock().unwrap();
    match map.remove(&id) {
        Some(sender) => {
            println!("[frontend_bool_response] Sending value: {}", value);
            let _ = sender.send(value);
        }
        None => {
            println!("[frontend_bool_response] No sender found for ID: {}", id);
            println!("Current keys: {:?}", map.keys().collect::<Vec<_>>());
        }
    }
}

#[tauri::command]
pub fn frontend_ready(app: AppHandle, state: State<CallbackMap>) {
    tauri::async_runtime::spawn({
        let app = app.clone();
        let state = state.inner().clone(); // clone the Arc<Mutex<...>>

        async move {
            println!("Checking for updates");
            info!("Checking for updates");
            update(app, state).await.unwrap();
        }
    });
}

pub async fn emit_and_wait_async(app: &tauri::AppHandle, state: CallbackMap) -> bool {
    let (tx, rx) = oneshot::channel();
    let id = uuid::Uuid::new_v4().to_string();

    {
        let mut map = state.lock().unwrap();
        map.insert(id.clone(), tx);
        println!("[emit_and_wait_async] Inserted sender for ID: {}", id);
    }

    app.emit("ask-frontend-bool", id.clone()).unwrap();
    println!("[emit_and_wait_async] Event emitted");

    match rx.await {
        Ok(result) => {
            println!("[emit_and_wait_async] Got result: {}", result);
            result
        }
        Err(e) => {
            println!("[emit_and_wait_async] ERROR: oneshot recv failed: {:?}", e);
            println!(
                "[emit_and_wait_async] Possibly sender was dropped before frontend responded."
            );
            false
        }
    }
}
/*
pub fn emit_and_wait(app: &AppHandle, state: CallbackMap) -> bool {
    let (tx, rx) = mpsc::channel();
    let id = Uuid::new_v4().to_string();

    // Store the sender with the unique ID
    state.lock().unwrap().insert(id.clone(), tx);
    // Emit event with the unique ID
    app.emit("ask-frontend-bool", id.clone()).unwrap();
    // Wait for response (blocking)
    let response: bool = rx.recv().unwrap();
    info!("Got response from frontend: {}", response);
    response
}
*/

pub async fn update(app: tauri::AppHandle, state: CallbackMap) -> tauri_plugin_updater::Result<()> {
    if let Some(update) = app.updater()?.check().await? {
        let allow_update = emit_and_wait_async(&app, state.clone()).await;
        println!("Confirmed from frontend: {}", allow_update); // <== make sure this prints!

        if allow_update {
            update
                .download_and_install(
                    |chunk, total| println!("Progress: {chunk}/{:?}", total),
                    || println!("Update complete"),
                )
                .await?;

            app.restart();
        } else {
            println!("User canceled update.");
        }
    } else {
        //let confirmed = emit_and_wait_async(&app, state.clone()).await;
        println!("Up-to-date");
    }

    Ok(())
}
