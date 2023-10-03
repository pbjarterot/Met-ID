use std::path::{Path, PathBuf};
use std::fs;
use log::info;
use crate::get_executable_path;



fn check_or_create(path: &PathBuf) {
  if !path.exists() {
    fs::create_dir_all(&path).expect("Failed to create AppData directory");
  }
}


fn copy_dir_to(src: &Path, dest: &Path) -> std::io::Result<()> {
  if !dest.exists() {
      fs::create_dir_all(&dest)?;
  }

  for entry in fs::read_dir(src)? {
      let entry = entry?;
      let entry_path = entry.path();
      let dest_child = dest.join(entry.file_name());

      if entry_path.is_dir() {
          copy_dir_to(&entry_path, &dest_child)?;
      } else {
          fs::copy(&entry_path, &dest_child)?;
      }
  }
  Ok(())
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

pub fn ensure_python_in_appdata(app: &tauri::App, resource_name: &str) -> std::path::PathBuf {
  // 1. Resolve the AppData path
  let mut app_data_path: PathBuf = app.path_resolver().app_data_dir().expect("Failed to resolve AppData path");
  app_data_path.push("met-id");
  check_or_create(&app_data_path);
  app_data_path.push("bin");

  check_or_create(&app_data_path);
  app_data_path.push(resource_name);

  info!("app_data_path for python is {:?}", app_data_path);


  let resource_path: PathBuf = get_executable_path().parent().unwrap().to_path_buf().join("bin");
  info!("resource_path for python is {:?}", resource_path);
  
  // Perform the actual copy
  fs::create_dir_all(app_data_path.parent().unwrap()).expect("Failed to create directories");
  info!("copying in: {:?}\n{:?}\n\n", resource_path, app_data_path);
  copy_dir_to(&resource_path, &app_data_path).expect("Failed to copy database");
  app_data_path
}