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


pub fn mass_matcher(input_masses: Vec<f64>, db_masses: &[f64], db_names: Vec<String>, db_mnames: Vec<String>, db_accessions: Vec<String>, db_smiles: Vec<String>) -> Vec<Vec<HashMap<String, String>>>{
    //input:
    //list of peaks to identify
    //database list of peaks, names and ms2 availability to identify against
    //mass error calibration

    let mut res: Vec<Vec<HashMap<String, String>>> = Vec::new();

    for mass in input_masses {
        let mut res_vec: Vec<HashMap<String, String>> = Vec::new();
        let low_mass: f64 = mass - 10.0*(mass/1E6);
        let high_mass: f64 = mass + 10.0*(mass/1E6);

        let low_idx: usize = bisect_left(db_masses, &low_mass);
        let high_idx: usize = bisect_left(db_masses, &high_mass);

        println!("{:?}, {:?}", low_mass, high_mass);
        println!("{:?}, {:?}", low_idx, high_idx);


        for index in low_idx..high_idx {
            let mut res_map: HashMap<String, String> = HashMap::new();
            res_map.insert("mz".to_string(), db_masses[index].to_string());
            res_map.insert("names".to_string(), db_names[index].to_string());
            res_map.insert("matrix".to_string(), db_mnames[index].to_string());
            res_map.insert("accession".to_string(), db_accessions[index].to_string());
            res_map.insert("smiles".to_string(), db_smiles[index].to_string());
            res_vec.push(res_map);
        }

        res.push(res_vec);
        
    };




    res
}