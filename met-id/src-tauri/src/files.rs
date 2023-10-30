use std::fs::File;
use std::io::prelude::*;
use serde::Deserialize;
use std::io;//::{self, BufRead};
use std::path::Path;
use std::fs;
use std::str::FromStr;
use csv::ReaderBuilder;

#[derive(Debug)]
enum Delim {
    Comma, 
    Tab,
    Semicolon,
}

impl FromStr for Delim {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "," => Ok(Delim::Comma),
            "\t" => Ok(Delim::Tab),
            ";" => Ok(Delim::Semicolon),
            _ => Err(()),
        }
    }
}

#[allow(dead_code)]
struct MS1File {
    filename: String,
    delimiter: Delim,
    mz: Vec<f64>,
    ion_mobility: Option<Vec<f64>>,
}

impl MS1File {
    fn parse_from_file(file_path: &Path) -> Result<Self, io::Error> {
        let file_content = fs::read_to_string(file_path)?;
        let lines: Vec<&str> = file_content.lines().collect();

        let filename = file_path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string();

        let delimiter_str = find_delimiter_line(&lines)?;
        let delimiter = parse_delimiter(&delimiter_str)?;

        let (mz, ion_mobility) = parse_data_points(&lines)?;

        Ok(Self {
            filename,
            delimiter,
            mz,
            ion_mobility,
        })
    }
}

fn find_delimiter_line<'a>(lines: &'a [&'a str]) -> Result<&'a str, io::Error> {
    lines.iter()
        .find(|&&line| !line.starts_with('#'))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "No delimiter line found")).copied()
}

fn parse_delimiter(delimiter_str: &str) -> Result<Delim, io::Error> {
    if delimiter_str.contains(';') {
        Ok(Delim::Semicolon)
    } else if delimiter_str.contains(',') {
        Ok(Delim::Comma)
    } else if delimiter_str.contains('\t') {
        Ok(Delim::Tab)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid delimiter: {}", delimiter_str),
        ))
    }
}

fn parse_data_points(lines: &[&str]) -> Result<(Vec<f64>, Option<Vec<f64>>), io::Error> {
    let data_lines: Vec<&str> = lines
        .iter()
        .filter(|&&line| !line.starts_with('#'))
        .skip(1)
        .cloned()
        .collect();

    let (mz, ion_mobility): (Vec<_>, Vec<_>) = data_lines
        .iter()
        .map(|line| {
            let values: Vec<&str> = line.split(';').map(str::trim).collect();
            let mz_value = values.get(0).and_then(|&v| v.parse().ok());
            let ion_mobility_value = values.get(1).and_then(|&v| v.parse().ok());
            (mz_value, ion_mobility_value)
        })
        .filter_map(|(mz_value, ion_mobility_value)| {
            mz_value.and_then(|mz: f64| ion_mobility_value.map(|ion_mobility: f64| (mz, ion_mobility)))
        })
        .unzip();

    let ion_mobility = if ion_mobility.is_empty() {
        None
    } else {
        Some(ion_mobility)
    };

    Ok((mz, ion_mobility))
}

#[derive(Deserialize)]
#[allow(non_snake_case)]
struct MSMSSpectra {
    name: String,
    matrix: String,
    SMILES: String,
    cid: String,
    fragment: String
}


#[tauri::command]
pub fn load_msms() -> Vec<Vec<std::string::String>> {
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


#[tauri::command]
pub fn parse_ms1_csv(path: String) -> Vec<f64>{
    let mut mzs : Vec<f64> = Vec::new();

    match MS1File::parse_from_file(&Path::new(&path)) {
        Ok(ms1_file) => {
            // Successfully parsed the file into an MS1File struct
            mzs = ms1_file.mz;
        }
        Err(err) => {
            // Failed to parse the file
            println!("Error: {}", err);
        }
    }
    mzs
    //Vec::new()
}

#[tauri::command]
pub fn read_ms1_ctrlv(content: String, delimiter: String, transpose: bool, _deletedrows: Vec<usize>) -> Vec<Vec<String>> {
    let lines: Vec<&str> = content.split("\n").collect();
    let mut res_rows: Vec<Vec<String>> = Vec::new();

    for line in lines {
        let res: Vec<String> = line.split(&delimiter).map(|s| s.to_string()).collect();
        res_rows.push(res);
    }

    if transpose {
        res_rows = transpose_vec(res_rows);
    }

    res_rows
}

#[tauri::command]
pub fn read_ms1_csv(path: String, delimiter: String, transpose: bool, _deletedrows: Vec<usize>) -> Vec<Vec<String>> {
    // Open the file
    let file: File = File::open(path).unwrap();

    // Create a CSV reader from the file
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
    let mut res_rows: Vec<Vec<String>> = Vec::new();

    for line in reader.records() {
        let rec: csv::StringRecord = line.unwrap();
        let record_str: &str = rec.get(0).unwrap();
        let res: Vec<String> = record_str.split(&delimiter).map(|s| s.to_string()).collect();
        res_rows.push(res);
    }

    // Transpose the dataset if `transpose` is true
    if transpose {
        res_rows = transpose_vec(res_rows);
    }
    /* 
    // Remove rows specified in `deleted_rows`
    let filtered_rows: Vec<Vec<String>> = res_rows
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !deleted_rows.contains(i))
        .map(|(_, row)| row)
        .collect();
    
    filtered_rows
    */

    res_rows
}



// Function to transpose a matrix
fn transpose_vec<T>(matrix: Vec<Vec<T>>) -> Vec<Vec<String>>
    where
        T: ToString,
    {
        // Find the length of the longest inner vector
        let max_len = matrix.iter().map(|row| row.len()).max().unwrap_or(0);

        // Create a new vector to store the transposed elements as strings
        let mut transposed: Vec<Vec<String>> = vec![Vec::new(); max_len];

        // Iterate through the rows of the original matrix
        for (row_index, row) in matrix.iter().enumerate() {
            // Iterate through the elements of each row and populate the corresponding columns in the transposed matrix
            for (col_index, elem) in row.iter().enumerate() {
                let value_str = elem.to_string();
                if transposed[col_index].len() <= row_index {
                    transposed[col_index].resize(row_index + 1, String::new());
                }
                transposed[col_index][row_index] = value_str;
            }
        }

        transposed
    }



#[tauri::command]
pub fn save_csv(path: String, csvcontent: String) {
    match fs::write(path, csvcontent){
        Ok(_) => println!("Success!"),
        Err(e) => println!("Failed to save csv: {:?}", e)
    };
}

/* 
#[tauri::command]
fn parse_text_file(file_path: &str) -> io::Result<Vec<f64>> {
    let path = Path::new(file_path);
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);

    let mut numbers = Vec::new();
    let mut is_first_row = true;

    for line in reader.lines() {
        if let Ok(line) = line {
            let values: Vec<f64> = line
                .split_whitespace()
                .map(|v| v.parse::<f64>())
                .filter_map(Result::ok)
                .collect();

            if is_first_row && values.len() > 1 {
                numbers.extend(values);
                break; // Stop processing rows
            } else if !is_first_row && values.len() == 1 {
                numbers.push(values[0]);
            }

            is_first_row = false;
        }
    }

    Ok(numbers)
}
*/


#[tauri::command]
pub fn read_mass_error_csv(path: String) -> Vec<Vec<String>> {
    // Open the file
    let file: File = File::open(path).unwrap();

    // Create a CSV reader from the file
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
    let mut res_rows: Vec<Vec<String>> = Vec::new();

    for (index, line) in reader.records().enumerate() {

        if index != 0 {
            let rec: csv::StringRecord = line.unwrap();

            let res: Vec<String> = vec![rec.get(0).unwrap().to_string(), rec.get(1).unwrap().to_string()];
            res_rows.push(res);
        }
        
    }

    res_rows
}
