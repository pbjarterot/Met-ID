fn main() {
    let _ = tonic_build::compile_protos("src/add_to_db/protos/streaming.proto");
    tauri_build::build()
}
