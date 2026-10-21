use tonic_prost_build::configure;

fn main() {
    configure()
        .out_dir("src/add_to_db")
        .compile_protos(
            &["src/add_to_db/protos/streaming.proto"],
            &["proto"])
            .expect("Failed to compile proto files");

    tauri_build::build()
}
