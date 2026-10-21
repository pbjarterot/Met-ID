use crate::get_app_handle;
use log::warn;
use std::fmt;
use tauri_plugin_shell::process::CommandEvent;
use tauri_plugin_shell::ShellExt;
use log::info;

// Custom error type
#[derive(Debug, Serialize)]
pub struct CommandError {
    message: String,
}



#[tauri::command]
pub fn sidecar_function(
    sidecar_name: String,
    sidecar_arguments: Vec<String>,
) -> Result<String, CommandError> {
    info!("Starting processing1...");
    let (mut rx, _child) = get_app_handle()
        .unwrap()
        .shell()
        .sidecar(sidecar_name)
        .expect("failed to create `my-sidecar` binary command")
        .args(sidecar_arguments)
        .spawn()
        .expect("Failed to spawn sidecar");

    let (tx, rx_data) = std::sync::mpsc::channel();
    info!("{:?}, {:?}", tx, rx_data);
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
    info!("Spawning sidecar...");

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



use log::{error};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc, Mutex,
    },
    time::Duration,
};
use tauri::{AppHandle, State};
use tauri_plugin_shell::{
    process::{CommandChild},
};
use tokio::sync::oneshot;


impl CommandError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for CommandError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message);
        Ok(())
    }
    
}

impl std::error::Error for CommandError {}

#[derive(Debug, Serialize)]
struct CcsRequest {
    id: u64,
    action: String,
    smiles: String,
    adduct: String,
}

#[derive(Debug, Deserialize)]
struct CcsResponse {
    id: u64,
    ok: bool,
    predicted_ccs: Option<f64>,
    mz: Option<f64>,
    error: Option<String>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CcsPredictionResult {
    pub predicted_ccs: f64,
    pub mz: f64,
}

type PendingResponse =
    oneshot::Sender<Result<CcsResponse, CommandError>>;



#[derive(Clone)]
pub struct CcsSidecarState {
    inner: Arc<CcsSidecarInner>,
}

impl Default for CcsSidecarState {
    fn default() -> Self {
        Self {
            inner: Arc::new(CcsSidecarInner {
                child: Mutex::new(None),
                pending: Mutex::new(HashMap::new()),
                next_id: AtomicU64::new(1),
            }),
        }
    }
}

use std::{
    io::{BufRead, BufReader, Write},
    path::Path,
    process::{Child, ChildStdin, Command, Stdio},
};

enum CcsChild {
    Tauri(CommandChild),
    Headless {
        // Keep the process handle alive while the API is running.
        _process: Child,
        stdin: ChildStdin,
    },
}

impl CcsChild {
    fn write_bytes(&mut self, bytes: &[u8]) -> Result<(), CommandError> {
        match self {
            Self::Tauri(child) => child.write(bytes).map_err(|error| {
                CommandError::new(format!(
                    "Could not write to CCS sidecar: {error}"
                ))
            }),
            Self::Headless { stdin, .. } => stdin.write_all(bytes).map_err(|error| {
                CommandError::new(format!(
                    "Could not write to CCS sidecar: {error}"
                ))
            }),
        }
    }
}

struct CcsSidecarInner {
    child: Mutex<Option<CcsChild>>,
    pending: Mutex<HashMap<u64, PendingResponse>>,
    next_id: AtomicU64,
}

pub fn start_ccs_sidecar(
    app: &AppHandle,
    state: &CcsSidecarState,
) -> Result<(), CommandError> {
    let command = app
        .shell()
        .sidecar("prediction_model")
        .map_err(|error| {
            CommandError::new(format!(
                "Could not create prediction sidecar: {error}"
            ))
        })?
        .arg("--sidecar-mode");

    let (mut receiver, child) = command.spawn().map_err(|error| {
        CommandError::new(format!(
            "Could not start prediction sidecar: {error}"
        ))
    })?;

    {
        let mut child_guard = state
            .inner
            .child
            .lock()
            .map_err(|_| CommandError::new("CCS child lock poisoned"))?;

        *child_guard = Some(CcsChild::Tauri(child));
    }

    let inner = Arc::clone(&state.inner);

    tauri::async_runtime::spawn(async move {
        while let Some(event) = receiver.recv().await {
            match event {
                CommandEvent::Stdout(bytes) => {
                    let line = String::from_utf8_lossy(&bytes);

                    let response =
                        match serde_json::from_str::<CcsResponse>(line.trim()) {
                            Ok(response) => response,
                            Err(error) => {
                                warn!(
                                    "Invalid CCS sidecar JSON: {error}; output={line}"
                                );
                                continue;
                            }
                        };

                    let sender = match inner.pending.lock() {
                        Ok(mut pending) => pending.remove(&response.id),
                        Err(_) => {
                            error!("CCS pending-response lock poisoned");
                            None
                        }
                    };

                    if let Some(sender) = sender {
                        let _ = sender.send(Ok(response));
                    } else {
                        warn!(
                            "No pending request for CCS response {}",
                            response.id
                        );
                    }
                }

                CommandEvent::Stderr(bytes) => {
                    let line = String::from_utf8_lossy(&bytes);
                    warn!("CCS sidecar: {}", line.trim());
                }

                CommandEvent::Error(message) => {
                    error!("CCS sidecar process error: {message}");
                }

                CommandEvent::Terminated(status) => {
                    let message = format!(
                        "CCS sidecar terminated with code {:?}",
                        status.code
                    );

                    error!("{message}");

                    if let Ok(mut child) = inner.child.lock() {
                        *child = None;
                    }

                    if let Ok(mut pending) = inner.pending.lock() {
                        let waiting_requests = std::mem::take(&mut *pending);

                        for (_, sender) in waiting_requests {
                            let _ = sender.send(Err(CommandError::new(
                                message.clone(),
                            )));
                        }
                    }

                    break;
                }

                _ => {}
            }
        }
    });

    Ok(())
}

#[tauri::command]
pub async fn sidecar_ccs(
    state: State<'_, Arc<CcsSidecarState>>,
    smiles: String,
    adduct: Option<String>,
) -> Result<CcsPredictionResult, CommandError> {
    //predict_with_sidecar(state.inner().as_ref(), smiles, adduct).await

    
    if smiles.trim().is_empty() {
        return Err(CommandError::new("SMILES cannot be empty"));
    }

    let id = state
        .inner
        .next_id
        .fetch_add(1, Ordering::Relaxed);

    let request = CcsRequest {
        id,
        action: "predict".to_string(),
        smiles,
        adduct: adduct.unwrap_or_else(|| "[M+H]+".to_string()),
    };

    let mut request_bytes =
        serde_json::to_vec(&request).map_err(|error| {
            CommandError::new(format!(
                "Could not serialize CCS request: {error}"
            ))
        })?;

    // One JSON object per line.
    request_bytes.push(b'\n');

    let (response_sender, response_receiver) =
        oneshot::channel();

    state
        .inner
        .pending
        .lock()
        .map_err(|_| {
            CommandError::new("CCS pending-response lock poisoned")
        })?
        .insert(id, response_sender);

    let write_result: Result<(), CommandError> = (|| {
        let mut child_guard = state
            .inner
            .child
            .lock()
            .map_err(|_| CommandError::new("CCS child lock poisoned"))?;

        let child = child_guard.as_mut().ok_or_else(|| {
            CommandError::new("CCS sidecar is not running")
        })?;

        child.write_bytes(&request_bytes)
    })();

    if let Err(error) = write_result {
        if let Ok(mut pending) = state.inner.pending.lock() {
            pending.remove(&id);
        }

        return Err(error);
    }

    let response = match tokio::time::timeout(
        Duration::from_secs(120),
        response_receiver,
    )
    .await
    {
        Ok(Ok(result)) => result?,

        Ok(Err(_)) => {
            return Err(CommandError::new(
                "CCS response channel was closed",
            ));
        }

        Err(_) => {
            if let Ok(mut pending) = state.inner.pending.lock() {
                pending.remove(&id);
            }

            return Err(CommandError::new(
                "CCS prediction timed out after 120 seconds",
            ));
        }
    };

    if !response.ok {
        return Err(CommandError::new(
            response
                .error
                .unwrap_or_else(|| "CCS prediction failed".to_string()),
        ));
    }

    let predicted_ccs = response.predicted_ccs.ok_or_else(|| {
        CommandError::new(
            "CCS sidecar returned no predicted_ccs value",
        )
    })?;
    let mz = response.mz.ok_or_else(|| {
        CommandError::new(
            "CCS sidecar returned no mz value"
        )
    })?;

    Ok(CcsPredictionResult {
        predicted_ccs,
        mz
    })
    
}

pub(crate) async fn predict_with_sidecar(
    state: &CcsSidecarState,
    smiles: String,
    adduct: Option<String>,
) -> Result<CcsPredictionResult, CommandError> {
    if smiles.trim().is_empty() {
        return Err(CommandError::new("SMILES cannot be empty"));
    }

    const ALLOWED_ADDUCTS: [&str; 4] =
        ["[M+H]+", "[M]+", "[M+Na]+", "[M+K]+"];

    let adduct = adduct
        .filter(|value| !value.trim().is_empty())
        .unwrap_or_else(|| "[M+H]+".to_string());

    if !ALLOWED_ADDUCTS.contains(&adduct.as_str()) {
        return Err(CommandError::new(format!(
            "Invalid adduct '{adduct}'. Allowed adducts: {}",
            ALLOWED_ADDUCTS.join(", ")
        )));
    }

    let id = state
        .inner
        .next_id
        .fetch_add(1, Ordering::Relaxed);

    let request = CcsRequest {
        id,
        action: "predict".to_string(),
        smiles,
        adduct: adduct,
    };

    let mut request_bytes =
        serde_json::to_vec(&request).map_err(|error| {
            CommandError::new(format!(
                "Could not serialize CCS request: {error}"
            ))
        })?;

    // One JSON object per line.
    request_bytes.push(b'\n');

    let (response_sender, response_receiver) =
        oneshot::channel();

    state
        .inner
        .pending
        .lock()
        .map_err(|_| {
            CommandError::new("CCS pending-response lock poisoned")
        })?
        .insert(id, response_sender);

    let write_result: Result<(), CommandError> = (|| {
        let mut child_guard = state
            .inner
            .child
            .lock()
            .map_err(|_| CommandError::new("CCS child lock poisoned"))?;

        let child = child_guard.as_mut().ok_or_else(|| {
            CommandError::new("CCS sidecar is not running")
        })?;

        child.write_bytes(&request_bytes)
    })();

    if let Err(error) = write_result {
        if let Ok(mut pending) = state.inner.pending.lock() {
            pending.remove(&id);
        }

        return Err(error);
    }

    let response = match tokio::time::timeout(
        Duration::from_secs(120),
        response_receiver,
    )
    .await
    {
        Ok(Ok(result)) => result?,

        Ok(Err(_)) => {
            return Err(CommandError::new(
                "CCS response channel was closed",
            ));
        }

        Err(_) => {
            if let Ok(mut pending) = state.inner.pending.lock() {
                pending.remove(&id);
            }

            return Err(CommandError::new(
                "CCS prediction timed out after 120 seconds",
            ));
        }
    };

    if !response.ok {
        return Err(CommandError::new(
            response
                .error
                .unwrap_or_else(|| "CCS prediction failed".to_string()),
        ));
    }

    let predicted_ccs = response.predicted_ccs.ok_or_else(|| {
        CommandError::new(
            "CCS sidecar returned no predicted_ccs value",
        )
    })?;
    let mz = response.mz.ok_or_else(|| {
        CommandError::new(
            "CCS sidecar returned no mz value"
        )
    })?;

    Ok(CcsPredictionResult {
        predicted_ccs,
        mz
    })
}

#[cfg(windows)]
use std::os::windows::process::CommandExt;

pub fn start_ccs_sidecar_process(
    state: &CcsSidecarState,
    executable_path: impl AsRef<Path>,
) -> Result<(), CommandError> {
    let mut command = Command::new(executable_path.as_ref());

    command
        .arg("--sidecar-mode")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null());

    // Prevent Windows from opening a console window for the sidecar.
    #[cfg(windows)]
    command.creation_flags(0x0800_0000); // CREATE_NO_WINDOW

    let mut process = command.spawn().map_err(|error| {
        CommandError::new(format!(
            "Could not start prediction sidecar: {error}"
        ))
    })?;

    let stdin = process.stdin.take().ok_or_else(|| {
        CommandError::new("Could not open CCS sidecar stdin")
    })?;

    let stdout = process.stdout.take().ok_or_else(|| {
        CommandError::new("Could not open CCS sidecar stdout")
    })?;

    {
        let mut child_guard = state
            .inner
            .child
            .lock()
            .map_err(|_| CommandError::new("CCS child lock poisoned"))?;

        *child_guard = Some(CcsChild::Headless {
            _process: process,
            stdin,
        });
    }

    let inner = Arc::clone(&state.inner);

    std::thread::spawn(move || {
        let reader = BufReader::new(stdout);

        for line in reader.lines() {
            match line {
                Ok(line) => {
                    let response =
                        match serde_json::from_str::<CcsResponse>(&line) {
                            Ok(response) => response,
                            Err(error) => {
                                warn!(
                                    "Invalid CCS sidecar JSON: {error}; output={line}"
                                );
                                continue;
                            }
                        };

                    let sender = match inner.pending.lock() {
                        Ok(mut pending) => pending.remove(&response.id),
                        Err(_) => {
                            error!("CCS pending-response lock poisoned");
                            None
                        }
                    };

                    if let Some(sender) = sender {
                        let _ = sender.send(Ok(response));
                    } else {
                        warn!(
                            "No pending request for CCS response {}",
                            response.id
                        );
                    }
                }

                Err(error) => {
                    error!("Could not read CCS sidecar output: {error}");
                    break;
                }
            }
        }

        let waiting_requests = match inner.pending.lock() {
            Ok(mut pending) => std::mem::take(&mut *pending),
            Err(_) => {
                error!("CCS pending-response lock poisoned");
                HashMap::new()
            }
        };

        for (_, sender) in waiting_requests {
            let _ = sender.send(Err(CommandError::new(
                "CCS sidecar closed its output",
            )));
        }

        if let Ok(mut child) = inner.child.lock() {
            *child = None;
        }
    });

    Ok(())
}