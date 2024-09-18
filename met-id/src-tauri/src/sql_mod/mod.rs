mod sql_handler;
mod build_query;
mod query;
mod counter;
mod msms;
mod functional_groups;
mod tissues;
mod matrices;
mod mtx_dropdown;
mod user_tables;
pub mod table;
mod database_view;

use self::database_view::{ get_db_data, db_ids_and_names};
use self::sql_handler::sql_handler;
use self::counter::sql_counter;
use self::msms::{ get_msms, get_msms_spectra, get_name_from_identifier_msms, ms2_search_spectra, match_msms_to_ui, add_msms_to_db, show_user_msms_db, remove_row_from_msms_user_db };
use self::functional_groups::get_functional_groups;
use self::tissues::get_tissues;
use self::matrices::get_matrices;
use self::user_tables::{ update_user_metabolites, remove_row_from_user_metabolites, remove_row_from_user_matrices, update_user_matrices, update_user_fgs, remove_user_fgs};
use self::mtx_dropdown::matrix_dropdown;

use std::collections::HashMap;
use serde::{ Serialize, Deserialize };



#[derive(Deserialize, Debug)]
pub struct Args {
    metabolome: String,
    matrix: String,
    met_type: Vec<String>,
    adducts: Vec<String>
}

#[derive(Debug)]
pub struct MS1DbRow {
    mz: f64,
    name: String,
    mname: String,
    accession: String,
    smiles: String,
    formula: String,
    possible_derivs: i32
}

#[derive(Debug)]
pub struct MS2DbRow {
    name: String,
    identifier: String,
    adduct: String,
    cid: String,
    window: String,
    tof: String,
    mz: String,
    data: (Vec<f64>, Vec<i64>),
    matrix: String
}

#[derive(Debug)]
pub struct MS2DbRow2 {
    id: usize,
    name: String,
    identifier: String,
    adduct: String,
    cid: String,
    window: String,
    tof: String,
    mz: String,
    data: (Vec<f64>, Vec<i64>),
    matrix: String
}

#[derive(Debug)]
pub struct IDentifierMap {
    name: String,
    hmdb: String,
}

#[derive(Serialize, Deserialize)]
pub struct SpectrumPoint {
    x: f64,
    y: i64,
}




#[tauri::command]
pub fn sql_handler_tauri(met: String, mat: String, typ: Vec<String>, adducts: Vec<String>, _mass_error: String, masses: Vec<String>, mzwindow: String) -> Vec<Vec<HashMap<String, String>>> {
    sql_handler(met, mat, typ, adducts, _mass_error, masses, mzwindow)
}

#[tauri::command]
pub fn sql_counter_tauri(met: String, mat: String, typ: Vec<String>, adducts: Vec<String>) -> i64 {
    sql_counter(met, mat, typ, adducts)
}

#[tauri::command]
pub fn get_msms_tauri() -> Vec<Vec<String>>{
    get_msms()
}

#[tauri::command]
pub fn get_msms_spectra_tauri(identifier: String, adduct: String, cid: String) -> String {
    get_msms_spectra(identifier, adduct, cid)
}

#[tauri::command]
pub fn get_name_from_identifier_msms_tauri(identifier: String) -> String {
    get_name_from_identifier_msms(identifier)
}

#[tauri::command]
pub fn ms2_search_spectra_tauri(name: String, fragment: String, ms1mass: String, fragmentslider: String, _ms1massslider: String) -> Vec<Vec<String>> {
    ms2_search_spectra(name, fragment, ms1mass, fragmentslider, _ms1massslider)
}

#[tauri::command]
pub fn get_functional_groups_tauri() -> Vec<String> {
    get_functional_groups()
}

#[tauri::command]
pub fn get_tissues_tauri() -> Vec<String> {
    get_tissues()
}

#[tauri::command]
pub fn get_matrices_tauri() -> Vec<String> {
    get_matrices()
}

#[tauri::command]
pub fn match_msms_to_ui_tauri(binsize: f64, threshold: f64) -> (Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<f64>, Vec<String>) {
    match_msms_to_ui(binsize, threshold)
}

#[tauri::command]
pub fn add_msms_to_db_tauri(name: String, adduct: String, mz: String, cid: String, tof: String, mzwindow: String, identifier: String, path: String, matrix: String) -> String {
    add_msms_to_db(name, adduct, mz, cid, tof, mzwindow, identifier, path, matrix)
}

#[tauri::command]
pub fn show_user_msms_db_tauri() -> Vec<Vec<String>> {
    show_user_msms_db()
}

#[tauri::command]
pub fn remove_row_from_msms_user_db_tauri(rowid: usize) -> usize {
    remove_row_from_msms_user_db(rowid)
}

#[tauri::command]
pub fn update_user_metabolites_tauri() -> Vec<Vec<String>> {
    update_user_metabolites()
}

#[tauri::command]
pub fn remove_row_from_user_metabolites_tauri(rowid: usize) -> usize {
    remove_row_from_user_metabolites(rowid)
}


#[tauri::command]
pub fn update_user_matrices_tauri() -> Vec<Vec<String>> {
    update_user_matrices()
}

#[tauri::command]
pub fn remove_row_from_user_matrices_tauri(rowid: usize) -> usize {
    remove_row_from_user_matrices(rowid)
}


#[tauri::command]
pub fn update_user_fgs_tauri() -> Vec<Vec<String>> {
    update_user_fgs()
}

#[tauri::command]
pub fn remove_from_user_fgs_tauri(rowid: usize, toremove: &str) -> usize {
    remove_user_fgs(rowid, toremove)
}

#[tauri::command]
pub fn matrix_dropdown_tauri() -> HashMap<String, Vec<String>> {
    matrix_dropdown()
}
#[tauri::command]
pub fn db_data_tauri(index: usize) -> (String, String, String, String, HashMap<String, HashMap<String, f64>>) {
    get_db_data(index)
}
#[tauri::command]
pub fn db_ids_and_names_tauri(inputvalue: String) -> Vec<(String, usize, i64)> {
    db_ids_and_names(inputvalue)
}

























