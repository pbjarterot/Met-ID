use std::collections::HashMap;

fn bisect_left<T: PartialOrd>(arr: &[T], x: &T) -> usize {
    let mut lo: usize = 0;
    let mut hi: usize = arr.len();
    while lo < hi {
        let mid: usize = (lo + hi) / 2;
        if arr[mid] < *x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}
#[derive(Debug, PartialEq)]
pub struct EquationCoefficients {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

pub fn parse_equation(input: &str) -> Option<EquationCoefficients> {
    let values: Vec<&str> = input.split(',').map(|s| s.trim()).collect();

    if input.is_empty() {
        return Some(EquationCoefficients { a: 0.0, b: 0.0, c: 0.0})
    }
    
    match values.len() {
        1 => {
            if let Ok(constant) = values[0].parse::<f64>() {
                Some(EquationCoefficients { a: 0.0, b: 0.0, c: constant })
            } else {
                None
            }
        },
        2 => {
            if let (Ok(a), Ok(b)) = (values[0].parse::<f64>(), values[1].parse::<f64>()) {
                Some(EquationCoefficients { a, b, c: 0.0 })
            } else {
                None
            }
        },
        3 => {
            if let (Ok(a), Ok(b), Ok(c)) = (values[0].parse::<f64>(), values[1].parse::<f64>(), values[2].parse::<f64>()) {
                Some(EquationCoefficients { a, b, c })
            } else {
                None
            }
        },
        _ => None,
    }
}

pub fn calculate_y(coefficients: &EquationCoefficients, x: f64) -> Option<f64> {
    let EquationCoefficients { a, b, c } = coefficients;
    let ppm_error: f64 = a * x * x + b * x + c;
    Some(x - (ppm_error*(x/1_000_000.0)))
}

fn msms_exist(msms_db: &Vec<String>, query: String) -> String {
    if msms_db.iter().any(|s| s == &query) {
        return "available".to_string()
    } else {
        return "not-available".to_string()
    }
}



fn update_hashmaps(data: &mut Vec<Vec<HashMap<String, String>>>) {
    let mut occurrences: HashMap<(String, String), i32> = HashMap::new(); // Map of (foo, bar) to its count
    let mut foo_counts: HashMap<String, i32> = HashMap::new(); // Map of foo to its total count

    // Count occurrences of each (foo, bar) pair and total foo counts
    for vec in data.iter() {
        for hashmap in vec.iter() {
            if let (Some(foo), Some(bar)) = (hashmap.get("names"), hashmap.get("matrix")) {
                *occurrences.entry((foo.clone(), bar.clone())).or_insert(0) += 1;
                *foo_counts.entry(foo.clone()).or_insert(0) += 1;
            }
        }
    }

    // Update the "coverage" value for each hashmap
    for vec in data.iter_mut() {
        for hashmap in vec.iter_mut() {
            if let Some(foo) = hashmap.get("names").map(String::as_str) {
                if let Some(total_foo_count) = foo_counts.get(foo) {
                    let mut count = *total_foo_count;
                    if let Some(bar) = hashmap.get("matrix").map(String::as_str) {
                        if let Some(same_pair_count) = occurrences.get(&(foo.to_string(), bar.to_string())) {
                            count -= same_pair_count-1;
                        }
                    }
                    hashmap.insert("coverage".to_string(), count.to_string());
                }
            }
        }
    }
}



#[tauri::command]
pub fn calculate_adjusted_mass(masses: Vec<String>, mass_error: String) -> Vec<HashMap<String, String>> {
    let mut input_masses: Vec<f64> = Vec::<f64>::new();
    let mut latent_variables: Vec<String> = Vec::<String>::new();
    for i in masses {
        match i.parse::<f64>() {
            Ok(num) => input_masses.push(num),
            Err(_e ) => {
                println!("Treating {:?} as a latent variable", i);
                latent_variables.push(i);
            },
        }
    }

    let mut res: Vec<HashMap<String, String>> = Vec::new();
    let coefficients: Option<EquationCoefficients> = parse_equation(&mass_error);
    for mass in input_masses {
        let adj_mass: f64;
        if let Some(ref coeff) = coefficients {
            if let Some(result) = calculate_y(&coeff, mass) {
                adj_mass = result;
            } else {
                println!("invalid_input");
                adj_mass = mass;
            }
        } else {
            println!("Invalid input");
            adj_mass = mass;
        }

        let mut res_map: HashMap<String, String> = HashMap::new();
        res_map.insert("oMass".to_string(), mass.to_string());
        res_map.insert("aMass".to_string(), adj_mass.to_string());
        res.push(res_map);
    }
    res
}


pub fn mass_matcher(input_masses: Vec<f64>, db_masses: &[f64], db_names: Vec<String>, db_mnames: Vec<String>, db_accessions: Vec<String>, db_smiles: Vec<String>, db_formulas: Vec<String>, db_pos_derivs: Vec<i32>, mass_error: String, mass_window: String, msms_ids: &Vec<String>) -> Vec<Vec<HashMap<String, String>>>{
    
    let mut res: Vec<Vec<HashMap<String, String>>> = Vec::new();
    let window: f64 = mass_window.parse::<f64>().unwrap();
    let coefficients: Option<EquationCoefficients> = parse_equation(&mass_error);

    for mass in input_masses {
        let mut res_vec: Vec<HashMap<String, String>> = Vec::new();
        let adj_mass: f64;
        if let Some(ref coeff) = coefficients {
            if let Some(result) = calculate_y(&coeff, mass) {
                adj_mass = result;
            } else {
                println!("invalid_input");
                adj_mass = mass;
            }
        } else {
            println!("Invalid input");
            adj_mass = mass;
        }

        
        let low_mass: f64 = adj_mass - (window/2.0)*(adj_mass/1E6);
        let high_mass: f64 = adj_mass + (window/2.0)*(adj_mass/1E6);

        let low_idx: usize = bisect_left(db_masses, &low_mass);
        let high_idx: usize = bisect_left(db_masses, &high_mass);



        for index in low_idx..high_idx {
            let db_name: String = if db_names[index].to_string() == "".to_string() {
                "Unknown Name".to_string()
            } else {
                db_names[index].to_string()
            };
            let mut res_map: HashMap<String, String> = HashMap::new();
            res_map.insert("oMass".to_string(), mass.to_string());
            res_map.insert("aMass".to_string(), adj_mass.to_string());
            res_map.insert("tMass".to_string(), db_masses[index].to_string());
            res_map.insert("dMass".to_string(), (adj_mass - db_masses[index]).to_string());
            res_map.insert("dPPM".to_string(), (((adj_mass - db_masses[index])/db_masses[index])*1E6).to_string());
            res_map.insert("names".to_string(), db_name);
            res_map.insert("matrix".to_string(), " [".to_string() + &db_mnames[index] + "]");
            res_map.insert("accession".to_string(), db_accessions[index].to_string());
            res_map.insert("smiles".to_string(), db_smiles[index].to_string());
            res_map.insert("msms".to_string(), msms_exist(msms_ids, db_accessions[index].to_string()));
            res_map.insert("formula".to_string(), db_formulas[index].to_string());
            res_map.insert("possible_derivs".to_string(), db_pos_derivs[index].to_string());
            res_vec.push(res_map);
        }

        if res_vec.is_empty() {
            let mut res_map: HashMap<String, String> = HashMap::new();
            res_map.insert("oMass".to_string(), mass.to_string());
            res_map.insert("aMass".to_string(), adj_mass.to_string());
            res_map.insert("tMass".to_string(), "".to_string());
            res_map.insert("dMass".to_string(), "".to_string());
            res_map.insert("dPPM".to_string(), "".to_string());
            res_map.insert("names".to_string(), "".to_string());
            res_map.insert("matrix".to_string(), "".to_string());
            res_map.insert("accession".to_string(), "".to_string());
            res_map.insert("smiles".to_string(), "".to_string());
            res_map.insert("msms".to_string(), "".to_string());
            res_map.insert("formula".to_string(), "".to_string());
            res_map.insert("possible_derivs".to_string(), "".to_string());
            res_vec.push(res_map);
        }
        res.push(res_vec);
        
    };

    update_hashmaps(&mut res);

    res
}