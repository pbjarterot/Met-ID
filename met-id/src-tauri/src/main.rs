// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod files;
mod sql;
mod mass_match;

fn main() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![files::load_msms, 
                                                 files::parse_ms1_csv, 
                                                 sql::sql_handler,
                                                 ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
