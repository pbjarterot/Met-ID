use std::fs::File;
use serde::Deserialize;
use std::io;
use std::path::Path;
use std::fs;
use std::str::FromStr;
use csv::ReaderBuilder;
use std::error::Error;
use imzml::mzml::{mzml::MzML, MzMLReader};
use std::io::{Read, Seek, SeekFrom};
use std::mem::size_of;
use imzml::{Compression, BinaryDataType};
use flate2::read::ZlibDecoder; 
use base64::{Engine as _, engine::general_purpose};
use byteorder::{ByteOrder, LittleEndian};
use std::io::Write;
//use crate::sidecar::sidecar_function;


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

/* 
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
*/

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


fn process_csv_data(reader: &mut csv::Reader<File> , delimiter: String) -> Result<Vec<Vec<String>>, Box<dyn Error>> {
    let mut res_rows = Vec::new();

    for result in reader.records() {
        match result {
            Ok(line) => {
                let rec: csv::StringRecord = line;
                if let Some(record_str) = rec.get(0) {
                    let res: Vec<String> = record_str.split(&delimiter.to_string()).map(|s| s.to_string()).collect();
                    res_rows.push(res);
                } else {
                    eprintln!("Warning: Empty record found.");
                }
            }
            Err(err) => {
                eprintln!("Warning: Error reading CSV line - Skipping: {}", err);
            }
        }
    }

    Ok(res_rows)
}


#[tauri::command]
pub fn read_ms1_csv(path: String, delimiter: String, transpose: bool, _deletedrows: Vec<usize>) -> Vec<Vec<String>> {
    // Open the file
    let file: File = File::open(path).unwrap();

    // Create a CSV reader from the file
    let mut reader: csv::Reader<File> = ReaderBuilder::new().has_headers(false).from_reader(file);
    //let res_rows: Vec<Vec<String>> = Vec::new();

    let mut res_rows: Vec<Vec<String>> = match process_csv_data(&mut reader, delimiter) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error: {}", e);
            return Vec::new() // or handle the error in some other way
        }
    };

    // Transpose the dataset if `transpose` is true
    if transpose {
        res_rows = transpose_vec(res_rows);
    }

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


fn read_characters(file_path: &str, start: u64, end: u64) -> io::Result<String> {
    let mut file = File::open(file_path)?;

    // Seek to the start position
    file.seek(SeekFrom::Start(start))?;

    // Calculate the number of bytes to read
    let length = end - start;

    // Read the specified number of bytes
    let mut buffer = vec![0; length as usize];
    file.read_exact(&mut buffer)?;

    // Convert the bytes to a String
    let content = String::from_utf8(buffer).expect("Invalid UTF-8 sequence");

    Ok(content)
}


fn parse_file(file_path: &str, offset: u64, length:u64, compression: Compression) -> io::Result<Vec<u8>> {
    let input_string: String = read_characters(file_path, offset, offset+length).unwrap();
    let bytes: Vec<u8> = general_purpose::STANDARD.decode(input_string).unwrap();

    let buffer = match compression {
        Compression::Zlib => {
            let mut decoder = ZlibDecoder::new(&bytes[..]);
            let mut decompressed_buffer = Vec::new();
            decoder.read_to_end(&mut decompressed_buffer)?;
            decompressed_buffer
        },
        Compression::None => {
            bytes.to_vec()
        },
        Compression::Undefined => {
            bytes.to_vec()
        },
    };

    Ok(buffer)
}


#[derive(Debug)]
#[derive(Clone)]
enum TypedData {
    F64Data(Vec<f64>),
    F32Data(Vec<f32>),
    U32Data(Vec<u32>),
}

fn convert_to_type(buffer: Vec<u8>, data_type: BinaryDataType) -> Result<TypedData, &'static str> {
    match data_type {
        BinaryDataType::Float64 => {
            if buffer.len() % size_of::<f64>() != 0 {
                return Err("Buffer length is not aligned with f64 size");
            }
            let f64_buffer: Vec<f64> = buffer.chunks_exact(size_of::<f64>())
                .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Ok(TypedData::F64Data(f64_buffer))
        },
        BinaryDataType::Float32 => {
            if buffer.len() % size_of::<f32>() != 0 {
                return Err("Buffer length is not aligned with u32 size");
            }
            let f32_buffer: Vec<f32> = buffer.chunks_exact(size_of::<f32>())
                .map(|chunk| f32::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Ok(TypedData::F32Data(f32_buffer))
        },
        BinaryDataType::Undefined => {
            if buffer.len() % size_of::<u32>() != 0 {
                return Err("Buffer length is not aligned with u32 size");
            }
            let u32_buffer: Vec<u32> = buffer.chunks_exact(size_of::<u32>())
                .map(|chunk| u32::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Ok(TypedData::U32Data(u32_buffer))
        },
    }
}



fn get_mzml_spectra(data: &imzml::BinaryDataArray, pth: &str) ->TypedData {
    let offset: u64 = data.offset().unwrap();
    let encoded_length: u64 = data.encoded_length().unwrap();
    let compression: imzml::Compression = data.compression();
    let type_: imzml::BinaryDataType = data.binary_type();

    let buffer: Result<Vec<u8>, io::Error> = parse_file(pth, offset, encoded_length, compression);
    let typed_data: TypedData = convert_to_type(buffer.unwrap(), type_).unwrap();

    //extract_and_print_first_10(&typed_data);
    //extract_and_print_details(&typed_data)
    typed_data
}




fn top_ten_thousand<T: PartialOrd + Copy>(vec: Vec<T>) -> Vec<usize> {
    let mut indexed_vec: Vec<(usize, T)> = vec.into_iter().enumerate().collect();
    
    // Sort by value in descending order
    indexed_vec.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Truncate to keep only top 1,000 elements
    if indexed_vec.len() > 100 {
        indexed_vec.truncate(100);
    }

    // Extract indices
    let indices: Vec<usize> = indexed_vec.into_iter().map(|(index, _)| index).collect();

    indices
}



fn extract_top_10(data: &TypedData, top_ten_thou_indices: &Vec<usize>) -> Vec<f64>{
    match data {
        TypedData::F64Data(vec) => {
            let vec2: Vec<f64> = top_ten_thou_indices.iter().map(|&index| vec[index]).collect();
            vec2
        }
        TypedData::F32Data(vec) => {
            let vec2: Vec<f32> = top_ten_thou_indices.iter().map(|&index| vec[index]).collect();
            vec2.into_iter().map(|x| x as f64).collect()

        }
        TypedData::U32Data(vec) => {
            let vec2: Vec<u32> = top_ten_thou_indices.iter().map(|&index| vec[index]).collect();
            vec2.into_iter().map(|x| x as f64).collect()
            
        }
    }
}


#[tauri::command]
pub fn read_mzml_for_msms(path: String) -> String {

    let (mz_data, intensity_data) = read_mzml_for_msms2(path);

    // Ensure both mz_data and intensity_data have the same length
    assert_eq!(mz_data.len(), intensity_data.len());

    // Serialize the data
    let mut blob: Vec<u8> = Vec::new();

    // Interleave mz_data and intensity_data
    for (&mz_value, &intensity_value) in mz_data.iter().zip(intensity_data.iter()) {
        // Serialize mz_value as f64
        let mut mz_buf: [u8; 8] = [0u8; 8];
        LittleEndian::write_f64(&mut mz_buf, mz_value);
        blob.write_all(&mz_buf).unwrap();

        // Serialize intensity_value as i64
        let mut intensity_buf: [u8; 8] = [0u8; 8];
        LittleEndian::write_i64(&mut intensity_buf, intensity_value as i64);
        blob.write_all(&intensity_buf).unwrap();
    }

    crate::add_to_db::add_to_db_functions::fill_user_msms(blob);

    String::from("Done")
}

pub fn read_mzml_for_msms_to_add_to_db(path: String) -> Vec<u8> {

    let (mz_data, intensity_data) = read_mzml_for_msms2(path);

    // Ensure both mz_data and intensity_data have the same length
    assert_eq!(mz_data.len(), intensity_data.len());

    // Serialize the data
    let mut blob: Vec<u8> = Vec::new();

    // Interleave mz_data and intensity_data
    for (&mz_value, &intensity_value) in mz_data.iter().zip(intensity_data.iter()) {
        // Serialize mz_value as f64
        let mut mz_buf: [u8; 8] = [0u8; 8];
        LittleEndian::write_f64(&mut mz_buf, mz_value);
        blob.write_all(&mz_buf).unwrap();

        // Serialize intensity_value as i64
        let mut intensity_buf: [u8; 8] = [0u8; 8];
        LittleEndian::write_i64(&mut intensity_buf, intensity_value as i64);
        blob.write_all(&intensity_buf).unwrap();
    }

    blob
}



pub fn read_mzml_for_msms2(path: String) -> (Vec<f64>, Vec<f64>){
    let pth: String = path.clone();
    let parser: MzMLReader<io::BufReader<File>> = MzMLReader::from_path(path).unwrap();

    let mzml: MzML = parser.into();

    let spectrum_list: &imzml::SpectrumList = mzml.spectrum_list().unwrap();
    let spec2: &std::sync::Arc<imzml::Spectrum> = spectrum_list.spectrum(0).unwrap();

    let mz: &imzml::BinaryDataArray = spec2.mz_array().unwrap();
    let intensity: &imzml::BinaryDataArray = spec2.intensity_array().unwrap();

    let mz_data: TypedData = get_mzml_spectra(mz, &pth[..]);
    let intensity_data: TypedData = get_mzml_spectra(intensity, &pth[..]);
    let also_intensity_data = intensity_data.clone();

    let top_ten_thou_indices: Vec<usize> = match also_intensity_data {
        TypedData::F64Data(vec) => {
            top_ten_thousand(vec)
        }
        TypedData::F32Data(vec) => {
            top_ten_thousand(vec)
        }
        TypedData::U32Data(vec) => {
            top_ten_thousand(vec)
        }
    };

    let top_ten_thou_intensity: Vec<f64> = extract_top_10(&intensity_data, &top_ten_thou_indices);
    let top_ten_thou_mz: Vec<f64> = extract_top_10(&mz_data, &top_ten_thou_indices);

    (top_ten_thou_mz, top_ten_thou_intensity)
}