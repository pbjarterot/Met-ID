pub mod errors;
pub mod identify;
pub mod sql;
pub mod parsing;
pub mod imzml;
pub mod ccs_prediction;

use axum::{routing::get, routing::post, Router, Json};
use std::net::SocketAddr;
use std::sync::Arc;


use crate::api::ccs_prediction::{ApiState, predict_ccs_api};
use crate::api::identify::{identify_metabolites, process_items};
use crate::api::sql::{custom_sql, get_database, post_database, search_by_name};
use crate::api::imzml::read_imzml;
use crate::sidecar::CcsSidecarState;


pub async fn serve_api(port: u16, ccs_state: Arc<CcsSidecarState>) {

    let router = Router::new()
        .route("/", get(|| async {"Hello from Tauri API!"}))
        .route("/warmup", get(|| async { "OK" }))
        .route("/name", get(return_name))
        .route("/process", post(process_items))
        .route("/identify", post(identify_metabolites))
        .route("/get_database", get(get_database))
        .route("/post_database", post(post_database))
        .route("/cutom_sql", post(custom_sql))
        .route("/search_by_name", get(search_by_name))
        .route("/read_imzml", post(read_imzml))
        .route("/predict_ccs", post(predict_ccs_api))
        .with_state(ApiState {ccs: ccs_state });

    let listener = tokio::net::TcpListener::bind(SocketAddr::from(([127, 0, 0, 1], port))).await.unwrap();
    axum::serve(listener, router).await.unwrap();
}

#[tauri::command]
pub async fn start_server(
    port: u16,
    ccs_state: tauri::State<'_, Arc<CcsSidecarState>>,
) -> Result<(), String> {
    serve_api(port, Arc::clone(ccs_state.inner())).await;
    Ok(())
}

async fn return_name(Json(name): Json<String>) -> Json<String> {
    Json(format!("Hello {}!", name))
}