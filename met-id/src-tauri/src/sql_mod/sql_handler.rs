use std::collections::HashMap;
use std::time::Instant;
use crate::mass_match::mass_matcher;

use super::build_query::build_query;
use super::Args;
use super::msms::get_msms;
use super::query::sql_query;



pub fn sql_handler(met: String, mat: String, typ: Vec<String>, adducts: Vec<String>, _mass_error: String, masses: Vec<String>, mzwindow: String) -> Vec<Vec<HashMap<String, String>>> {
  let count: bool = false;

    let args: Args = Args {
        metabolome: met,
        matrix: mat,
        met_type: typ,
        adducts
    };

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

    let query_str: String = build_query(&args, min_peak - 1.0, max_peak + 1.0, count, "".to_string());
    let query_str2: String = build_query(&args, min_peak - 1.0, max_peak + 1.0, count, "user_".to_string());
    println!("{:?}", query_str);

    let mut db_input: (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<i32>) = sql_query(&query_str);
    let db_input2: (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<i32>) = sql_query(&query_str2);
    println!("{:?}", db_input);

    db_input = merge_sorted_tuples(db_input, db_input2);

    let msms_ids: Vec<String> = get_msms_ids();

    
    let db_input_array: Box<[f64]> = db_input.0.into_boxed_slice();
    let db_input_array_ref: &[f64] = &*db_input_array;


    let start = Instant::now();                                      
    let output: Vec<Vec<HashMap<String, String>>> = mass_matcher(input_masses, db_input_array_ref, db_input.1, db_input.2, db_input.3, db_input.4, db_input.5, db_input.6, _mass_error, mzwindow, &msms_ids);    
    let duration = start.elapsed();
    println!("Time elapsed in ms1 matcher is: {:?}", duration);
    
    output
}

fn get_msms_ids() -> Vec<String> {
  let ids: Vec<String> = get_msms()[1].clone();
  ids
}



fn merge_sorted_tuples(
    tuple1: (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<i32>),
    tuple2: (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<i32>),
) -> (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<i32>) {
    let mut result = (
        Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()
    );

    let mut iter1 = (
        tuple1.0.into_iter(),
        tuple1.1.into_iter(),
        tuple1.2.into_iter(),
        tuple1.3.into_iter(),
        tuple1.4.into_iter(),
        tuple1.5.into_iter(),
        tuple1.6.into_iter(),
    );

    let mut iter2 = (
        tuple2.0.into_iter(),
        tuple2.1.into_iter(),
        tuple2.2.into_iter(),
        tuple2.3.into_iter(),
        tuple2.4.into_iter(),
        tuple2.5.into_iter(),
        tuple2.6.into_iter(),
    );

    let mut next1 = (
        iter1.0.next(),
        iter1.1.next(),
        iter1.2.next(),
        iter1.3.next(),
        iter1.4.next(),
        iter1.5.next(),
        iter1.6.next(),
    );
    
    let mut next2 = (
        iter2.0.next(),
        iter2.1.next(),
        iter2.2.next(),
        iter2.3.next(),
        iter2.4.next(),
        iter2.5.next(),
        iter2.6.next(),
    );

    while let (Some(v1), Some(v2)) = (next1.0, next2.0) {
        if v1 <= v2 {
            result.0.push(v1);
            result.1.push(next1.1.take().unwrap());
            result.2.push(next1.2.take().unwrap());
            result.3.push(next1.3.take().unwrap());
            result.4.push(next1.4.take().unwrap());
            result.5.push(next1.5.take().unwrap());
            result.6.push(next1.6.take().unwrap());
            next1 = (
                iter1.0.next(),
                iter1.1.next(),
                iter1.2.next(),
                iter1.3.next(),
                iter1.4.next(),
                iter1.5.next(),
                iter1.6.next(),
            );
        } else {
            result.0.push(v2);
            result.1.push(next2.1.take().unwrap());
            result.2.push(next2.2.take().unwrap());
            result.3.push(next2.3.take().unwrap());
            result.4.push(next2.4.take().unwrap());
            result.5.push(next2.5.take().unwrap());
            result.6.push(next2.6.take().unwrap());
            next2 = (
                iter2.0.next(),
                iter2.1.next(),
                iter2.2.next(),
                iter2.3.next(),
                iter2.4.next(),
                iter2.5.next(),
                iter2.6.next(),
            );
        }
    }

    while let Some(v1) = next1.0 {
        result.0.push(v1);
        result.1.push(next1.1.take().unwrap());
        result.2.push(next1.2.take().unwrap());
        result.3.push(next1.3.take().unwrap());
        result.4.push(next1.4.take().unwrap());
        result.5.push(next1.5.take().unwrap());
        result.6.push(next1.6.take().unwrap());
        next1 = (
            iter1.0.next(),
            iter1.1.next(),
            iter1.2.next(),
            iter1.3.next(),
            iter1.4.next(),
            iter1.5.next(),
            iter1.6.next(),
        );
    }

    while let Some(v2) = next2.0 {
        result.0.push(v2);
        result.1.push(next2.1.take().unwrap());
        result.2.push(next2.2.take().unwrap());
        result.3.push(next2.3.take().unwrap());
        result.4.push(next2.4.take().unwrap());
        result.5.push(next2.5.take().unwrap());
        result.6.push(next2.6.take().unwrap());
        next2 = (
            iter2.0.next(),
            iter2.1.next(),
            iter2.2.next(),
            iter2.3.next(),
            iter2.4.next(),
            iter2.5.next(),
            iter2.6.next(),
        );
    }

    result
}


