use tonic_prost_build::configure;

fn main() {
    configure()
        .out_dir("src/add_to_db")
        .compile_protos(
            &["src/add_to_db/protos/streaming.proto"],
            &["proto"])
            .expect("Failed to compile proto files");

    let attributes = tauri_build::Attributes::new().windows_attributes(
        tauri_build::WindowsAttributes::new().static_vc_runtime(true),
    );

    tauri_build::try_build(attributes)
        .expect("Failed to build Tauri application");
}
