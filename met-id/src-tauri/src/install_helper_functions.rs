use log::info;
use r2d2::Pool;
use r2d2_sqlite::SqliteConnectionManager;
use std::fs;
use std::path::PathBuf;
use tauri::App;
use tauri::Manager;

use crate::get_app_handle;

pub fn ensure_database_in_appdata(_app: &App, resource_name: &str) -> std::path::PathBuf {
    // Resolve the appropriate App Data path for the OS
    let mut app_data_path: PathBuf = get_app_handle()
        .unwrap()
        .path()
        .app_data_dir()
        .expect("Failed to resolve AppData path");

    // Create directories if they don't exist
    let dirs_to_create: [&str; 2] = ["met-id", "DBs"];
    for dir in &dirs_to_create {
        app_data_path.push(dir);
        if !app_data_path.exists() {
            fs::create_dir_all(&app_data_path).expect("Failed to create directory");
        }
    }
    app_data_path.push(resource_name);

    // Check if the database exists in AppData, if not, copy from resources
    if !app_data_path.exists() {
        let resource_path: PathBuf = get_app_handle()
            .unwrap()
            .path()
            .resource_dir()
            .expect("Failed to resolve resource directory")
            .join("DBs")
            .join(resource_name);

        fs::create_dir_all(app_data_path.parent().unwrap()).expect("Failed to create directories");
        fs::copy(&resource_path, &app_data_path).expect("Failed to copy database");
    }
    app_data_path
}

/*
pub fn ensure_bin_in_appdata(app: &App, resource_name: &str) -> std::path::PathBuf {
    // Resolve the appropriate App Data path for the OS
    let mut app_data_path: PathBuf = app.path_resolver().app_data_dir()
    .expect("Failed to resolve AppData path");

    // Create directories if they don't exist
    let dirs_to_create: [&str; 2] = ["met-id", "bin"];
    for dir in &dirs_to_create {
    app_data_path.push(dir);
    if !app_data_path.exists() {
        fs::create_dir_all(&app_data_path).expect("Failed to create directory");
    }
    }
    app_data_path.push(resource_name);

    // Check if the database exists in AppData, if not, copy from resources
    if !app_data_path.exists() {
    let resource_path: PathBuf = app.path_resolver().resource_dir().expect("Failed to resolve resource directory")
        .join("pyinstaller")
        .join("dist")
        .join(resource_name);

    println!("{:?}\n{:?}", app_data_path, resource_path);
    fs::create_dir_all(app_data_path.parent().unwrap()).expect("Failed to create directories");
    fs::copy(&resource_path, &app_data_path).expect("Failed to copy database");
    }
    app_data_path
}
*/

pub fn create_pool_from_app_path(
    _app: &App,
    resource_path: PathBuf,
) -> Pool<SqliteConnectionManager> {
    info!("looking in: {:?}", resource_path);
    println!("looking in: {:?}", resource_path);
    let manager = SqliteConnectionManager::file(resource_path);
    Pool::new(manager).expect("Failed to create pool.")
}
