// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
use std::fs::File;
use std::io::prelude::*;
use serde::{Deserialize};
use csv;

// Import the necessary Tauri and Rust modules

#[derive(Deserialize)]
struct MSMSSpectra {
    name: String,
    matrix: String,
    SMILES: String,
    cid: String,
    fragment: String

}

#[tauri::command]
fn load_msms() -> Vec<Vec<std::string::String>> {
    //needed: name, adduct, matrix, parent mass, smiles, spectra

    let mut file: File = File::open("../databases/sample2.json").expect("Failed to open file");

    // Read the contents of the file into a string
    let mut contents: String = String::new();
    file.read_to_string(&mut contents).expect("Failed to read file");

    // Deserialize the JSON data into your struct
    let data: Vec<MSMSSpectra> = serde_json::from_str(&contents).expect("Failed to parse JSON");

    let names: Vec<String> = data.iter().map(|x: &MSMSSpectra| x.name.clone()).collect();
    let _matrices: Vec<String> = data.iter().map(|x: &MSMSSpectra| x.matrix.clone()).collect();
    let smiles: Vec<String> = data.iter().map(|x: &MSMSSpectra| x.SMILES.clone()).collect();
    let cid: Vec<String> = data.iter().map(|x: &MSMSSpectra| x.cid.clone()).collect();
    let fragment: Vec<String> = data.iter().map(|x: &MSMSSpectra| x.fragment.clone()).collect();

    let return_data: Vec<Vec<String>> = vec![names, smiles, cid, fragment];
    // Use the data as needed
    return_data

}

fn check_delimiter(line: &str) -> &str {
    if line.contains(';') && !line.contains(",") && !line.contains("\t") {
        ";"
    } else if !line.contains(';') && line.contains(",") && !line.contains("\t") {
        ","
    } else if !line.contains(';') && !line.contains(",") && line.contains("\t") {
        "\t"
    } else {
        ","
    }
    
}

fn get_mz_from_record(input: &str, delimiter: &str) -> Option<f64> {
    let parts: Vec<&str> = input.split(delimiter).collect();
    let first_part: String = parts[0].to_string();

    match first_part.parse::<f64>(){
        Ok(num) => Some(num),
        Err(_) => {
            None
        },
    }
}



#[tauri::command]
fn parse_ms1_csv(path: String) -> Vec<f64>{
    let file: File = File::open(path).unwrap();

    // Create a CSV reader from the file
    let mut reader: csv::Reader<File> = csv::Reader::from_reader(file);
    let mut mzs: Vec<f64> = Vec::new();
    let mut index: usize = 0;

    for line in reader.records() {
        let record: csv::StringRecord = line.unwrap();
        let record_str: &str = record.get(0).unwrap();

        //CSV files from scils have some # values at the top so we dont want to parse that
        if record_str.starts_with("#") {
            continue
        }

        //in the first row of data we check whether the CSV is ';' delimited, ',' delimited or 'tab' delimited.
        let mut delimiter: &str = ",";
        if index == 0 {
            delimiter = check_delimiter(record_str);
        } 

        //parsing the data
        let mz: Option<f64> = get_mz_from_record(record_str, delimiter);
        if mz != None {
            mzs.push(mz.unwrap());
        }
        index += 1;
    }
    mzs
}


// Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
#[tauri::command]
fn greet(name: &str) -> () {
    println!("Hello, {}! You've been greeted from Rust!", name)
}

fn main() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![greet, load_msms, parse_ms1_csv])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
