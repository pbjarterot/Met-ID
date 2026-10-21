use log::info;

use crate::mass_match::mass_matcher;
use crate::sql_mod::query::QueryResults;
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
    ccss: Vec<String>,
    mzwindow: f64,
    ion_mobility: bool
) -> Vec<Vec<HashMap<&'static str, String>>> {
    let count: bool = false;
    let start2: Instant = Instant::now();
    let args: Args = Args {
        metabolome: met,
        matrix: mat,
        met_type: typ,
        adducts,
    };

    let ccs_included: bool = if ccss[0] != "-1" && (masses.len() == ccss.len()) {
        true
    } else {
        false
    };

    let mut input_masses: Vec<f64> = Vec::<f64>::new();
    let mut input_ccss: Vec<f64> = Vec::<f64>::new();
    let mut latent_variables: Vec<String> = Vec::<String>::new();

    for (idx, i) in masses.iter().enumerate() {
        match i.parse::<f64>() {
            Ok(num) => {
                input_masses.push(num);
                if ccs_included {
                    match ccss[idx].parse::<f64>() {
                        Ok(n) => input_ccss.push(n),
                        Err(_e) => {
                            println!("{:?} is not a valid ccs value", ccss[idx]);
                            input_ccss.push(0.0);
                        }
                    };
                }
            }
            Err(_e) => {
                println!("Treating {:?} as a latent variable", i);
                latent_variables.push(i.to_string());
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

    let query_str: String = build_query(&args, min_peak - 1.0, max_peak + 1.0, count, ion_mobility);
    let db_input: QueryResults = sql_query(&query_str);

    let msms_ids: Vec<String> = get_msms_ids();

    let db_input_array: Box<[f64]> = db_input.db_masses.into_boxed_slice();
    let db_input_array_ref: &[f64] = &*db_input_array;

    let db_input_array_ccs: Box<[f64]> = db_input.db_ccs.into_boxed_slice();
    let db_input_array_ref_ccs: &[f64] = &*db_input_array_ccs;

    let start = Instant::now();

    let output: Vec<Vec<HashMap<&'static str, String>>> = mass_matcher(
        input_masses,
        input_ccss,
        db_input_array_ref,
        db_input_array_ref_ccs, // to be changed when the database has ccs
        db_input.db_names,
        db_input.db_mnames,
        db_input.db_accessions,
        db_input.db_smiles,
        db_input.db_formulas,
        db_input.db_pos_derivs,
        _mass_error,
        mzwindow,
        &msms_ids,
        ccs_included,
    );
    let duration = start.elapsed();
    info!("Time elapsed in ms1 matcher is: {:?}", duration);

    let duration2 = start2.elapsed();
    info!("Full time is: {:?}", duration2);
    output
}

fn get_msms_ids() -> Vec<String> {
    let ids: Vec<String> = get_msms()[1].clone();
    ids
}
