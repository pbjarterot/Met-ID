use crate::mass_match::mass_matcher;
use crate::sql_mod::table::check_if_table_exists;
use std::collections::HashMap;
use std::time::Instant;

use super::build_query::build_query;
use super::msms::get_msms;
use super::query::sql_query;
use super::Args;

pub fn sql_handler(
    met: String,
    mat: String,
    typ: Vec<String>,
    adducts: Vec<String>,
    _mass_error: String,
    masses: Vec<String>,
    mzwindow: String,
) -> Vec<Vec<HashMap<String, String>>> {
    let count: bool = false;
    let start2 = Instant::now();
    let args: Args = Args {
        metabolome: met,
        matrix: mat,
        met_type: typ,
        adducts,
    };

    let mut input_masses: Vec<f64> = Vec::<f64>::new();
    let mut latent_variables: Vec<String> = Vec::<String>::new();
    for i in masses {
        match i.parse::<f64>() {
            Ok(num) => input_masses.push(num),
            Err(_e) => {
                println!("Treating {:?} as a latent variable", i);
                latent_variables.push(i);
            }
        }
    }

    let mut min_peak: f64 = 5000.0;
    let mut max_peak: f64 = 0.0;

    for &num in &input_masses {
        if num < min_peak {
            min_peak = num;
        }
        if num > max_peak {
            max_peak = num;
        }
    }

    check_if_table_exists("metabolites", "user_metabolites").unwrap();
    check_if_table_exists("derivatized_by", "user_derivatized_by").unwrap();
    check_if_table_exists("endogeneity", "user_endogeneity").unwrap();
    check_if_table_exists("in_tissue", "user_in_tissue").unwrap();
    check_if_table_exists("db_accessions", "user_db_accessions").unwrap();
    check_if_table_exists("functional_groups", "user_functional_groups").unwrap();
    let query_str: String = build_query(&args, min_peak - 1.0, max_peak + 1.0, count);
    println!("{:?}", query_str);
    let db_input: (
        Vec<f64>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<i32>,
    ) = sql_query(&query_str);

    println!("masses: {:?}", db_input.0);
    let msms_ids: Vec<String> = get_msms_ids();

    let db_input_array: Box<[f64]> = db_input.0.into_boxed_slice();
    let db_input_array_ref: &[f64] = &*db_input_array;

    let start = Instant::now();
    let output: Vec<Vec<HashMap<String, String>>> = mass_matcher(
        input_masses,
        db_input_array_ref,
        db_input.1,
        db_input.2,
        db_input.3,
        db_input.4,
        db_input.5,
        db_input.6,
        _mass_error,
        mzwindow,
        &msms_ids,
    );
    let duration = start.elapsed();
    println!("Time elapsed in ms1 matcher is: {:?}", duration);

    let duration2 = start2.elapsed();
    println!("Full time is: {:?}", duration2);
    output
}

fn get_msms_ids() -> Vec<String> {
    let ids: Vec<String> = get_msms()[1].clone();
    ids
}
