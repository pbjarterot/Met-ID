use log::{info, warn};
use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::params;
use std::sync::Arc;
use streaming::streaming_service_server::{StreamingService, StreamingServiceServer};
use streaming::MessageBatch;
use tokio::sync::{mpsc, Mutex};
use tokio_stream::wrappers::ReceiverStream;
use tonic::{transport::Server, Response, Status};

use crate::add_to_db::functional_group::get_smiles_from_db;
use crate::database_setup::get_connection;
use crate::get_app_handle;
use crate::sidecar::sidecar_function3;

pub mod streaming {
    //tonic::include_proto!("./streaming");
    include!("streaming.rs");
}

#[derive(Debug)]
pub struct MyStreamingService {
    shutdown_tx: Arc<Mutex<Option<tokio::sync::oneshot::Sender<()>>>>, // For graceful shutdown
    table_name: String,
    column_name: String,
    progress_sender: std::sync::mpsc::Sender<f32>,
    smiles_table: String
}

#[tonic::async_trait]
impl StreamingService for MyStreamingService {
    type StreamBatchedMessagesStream = ReceiverStream<Result<MessageBatch, Status>>;

    async fn stream_batched_messages(
        &self,
        request: tonic::Request<tonic::Streaming<MessageBatch>>,
    ) -> Result<Response<Self::StreamBatchedMessagesStream>, Status> {
        let mut conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

        let mut stream = request.into_inner();
        let (tx, rx) = mpsc::channel(4);

        let shutdown_tx = Arc::clone(&self.shutdown_tx);

        let table_name: String = self.table_name.clone();
        let column_name: String = self.column_name.clone();
        let progress_sender: std::sync::mpsc::Sender<f32> = self.progress_sender.clone();
        let smiles_table = self.smiles_table.clone();

        tokio::spawn(async move {
            // Step 1: Send each message individually
            let mut sent_count = 0;
            let mut received_count = 0;

            let messages: Vec<Vec<String>> = get_smiles_from_db(smiles_table);

            for message in messages {
                let message_batch = MessageBatch { messages: message };

                println!("Sending message to client: {:?}", sent_count);
                info!("Sending message to client: {:?}", sent_count);
                sent_count = sent_count + 1;
                progress_sender
                    .send(((sent_count as f32) / 220.0) * 50.0)
                    .unwrap();
                if let Err(e) = tx.send(Ok(message_batch)).await {
                    eprintln!("Failed to send message: {}", e);
                    warn!("Failed to send message: {}", e);
                    return;
                }

                //tokio::time::sleep(tokio::time::Duration::from_secs(0.01)).await; // Optional delay
            }

            // Step 2: Server sends stop signal
            let stop_message = MessageBatch {
                messages: vec!["stop".to_string()],
            };

            println!("Sending stop message to client: {:?}", stop_message);
            info!("Sending stop message to client: {:?}", stop_message);
            if let Err(e) = tx.send(Ok(stop_message)).await {
                eprintln!("Failed to send stop message: {}", e);
                warn!("Failed to send stop message: {}", e);
                return;
            }

            // Step 3: Wait for client's stop response
            while let Ok(Some(batch)) = stream.message().await {
                //println!("Received batch from client: {:?}", batch.messages);

                if batch.messages.contains(&"stop".to_string()) {
                    println!("Received stop confirmation from client. Shutting down...");
                    info!("Received stop confirmation from client. Shutting down...");
                    let mut shutdown_tx = shutdown_tx.lock().await;
                    if let Some(tx) = shutdown_tx.take() {
                        let _ = tx.send(()); // Trigger the shutdown signal
                    }
                    break;
                } else {
                    let transaction = conn.transaction().unwrap();

                    for (idx, message) in batch.messages[0].split(", ").enumerate() {
                        let new_idx: usize = (received_count * 1_000) + idx + 1;

                        transaction
                            .execute(
                                &format!("UPDATE '{}' SET '{}' = ?1 WHERE rowid = ?2", table_name, column_name),
                                params![message, new_idx],
                            )
                            .unwrap();

                        if idx == 0 {
                            println!("message: {:?}, idx: {:?}", message, new_idx);
                            info!("message: {:?}, idx: {:?}", message, new_idx);
                        }
                    }
                    transaction.commit().unwrap();
                }
                progress_sender
                    .send((((received_count as f32) / 220.0) * 50.0) + 50.0)
                    .unwrap();
                received_count = received_count + 1;
            }
            println!("sent: {:?}, received: {:?}", sent_count, received_count);
            println!("Stream handling completed.");
            info!("sent: {:?}, received: {:?}", sent_count, received_count);
            info!("Stream handling completed.");
        });
        Ok(Response::new(ReceiverStream::new(rx)))
    }
}

#[tokio::main]
pub async fn functional_group_server(
    progress_sender: std::sync::mpsc::Sender<f32>,
    smarts: String,
    table_name: String,
    column_name: String,
    smiles_table: String
) -> Result<(), Box<dyn std::error::Error>> {
    let (shutdown_tx, shutdown_rx) = tokio::sync::oneshot::channel();

    let addr = "[::1]:50051".parse()?;
    let streaming_service = MyStreamingService {
        shutdown_tx: Arc::new(Mutex::new(Some(shutdown_tx))),
        table_name: table_name.clone(),
        column_name: column_name.clone(),
        progress_sender: progress_sender.clone(),
        smiles_table: smiles_table.clone()
    };
    println!(
        "table_name: {:?}, column_name: {:?}, smarts: {:?}",
        table_name, column_name, smarts
    );
    info!(
        "Table name: {:?}, column name: {:?}",
        table_name, column_name
    );

    println!("Server listening on {}", addr);
    info!("Server listening on {}", addr);
    let sidecar_arguments: Vec<String> = vec![smarts];

    sidecar_function3(
        get_app_handle().unwrap(),
        progress_sender.clone(),
        sidecar_arguments,
    )
    .unwrap();

    // Graceful shutdown
    Server::builder()
        .add_service(StreamingServiceServer::new(streaming_service))
        .serve_with_shutdown(addr, async {
            shutdown_rx.await.ok();
        })
        .await?;

    println!("Server has shut down gracefully.");
    info!("Server has shut down gracefully.");
    Ok(())
}
