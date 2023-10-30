use std::path::PathBuf;
use std::fs;
use log::info;
use r2d2::Pool;
use r2d2_sqlite::SqliteConnectionManager;
use tauri::App;
use std::env;


fn check_or_create(path: &PathBuf) {
  if !path.exists() {
    fs::create_dir_all(&path).expect("Failed to create AppData directory");
  }
}



pub fn ensure_database_in_appdata(app: &tauri::App, resource_name: &str) -> std::path::PathBuf {
  // 1. Resolve the AppData path
  let mut app_data_path: PathBuf = app.path_resolver().app_data_dir().expect("Failed to resolve AppData path");
  check_or_create(&app_data_path);
  app_data_path.push("met-id");
  check_or_create(&app_data_path);
  app_data_path.push("DBs");
  check_or_create(&app_data_path);
  app_data_path.push(resource_name);

  // 2. Check if the database exists in AppData
  if !app_data_path.exists() {
      let resource_path: PathBuf = get_executable_path().parent().unwrap().to_path_buf().join("DBs").join(resource_name);


      // Perform the actual copy
      fs::create_dir_all(app_data_path.parent().unwrap()).expect("Failed to create directories");
      fs::copy(&resource_path, &app_data_path).expect("Failed to copy database");
  }
  app_data_path
}



pub fn create_pool_from_app_path(_app: &App, resource_path: PathBuf) -> Pool<SqliteConnectionManager> {
  info!("looking in: {:?}", resource_path);
  let manager = SqliteConnectionManager::file(resource_path);
  Pool::new(manager).expect("Failed to create pool.")
}

pub fn get_executable_path() -> std::path::PathBuf {
  env::current_exe().expect("Failed to get current executable path")
}