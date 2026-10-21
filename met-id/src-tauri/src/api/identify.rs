use serde::{Deserialize, Serialize};

use crate::{
    api::{
        errors::{APIResult, MetIDParseError},
        parsing::{
            strnum_to_vec_string, to_num, to_vec_string, NumInput, StringInput,
        },
    },
    sql_mod::sql_handler::sql_handler,
};

use std::collections::HashMap;
use axum::Json;
use super::errors::JsonParseError;

#[derive(Debug, Deserialize)]
pub struct InputData(Vec<f32>);

#[derive(Debug, Serialize)]
pub struct OutputData(Vec<f32>);

#[derive(Debug, Deserialize)]
pub struct IdentificationInputData {
    met: String,
    mat: String,
    typ: StringInput,
    adducts: StringInput,
    mass_error: NumInput,
    masses: NumInput,
    ccss: Option<NumInput>,
    mzwindow: NumInput,
}

#[derive(Debug, Serialize)]
pub struct IdentificationOutputData(HashMap<&'static str, Vec<String>>);

pub async fn identify_metabolites(
    Json(data): Json<IdentificationInputData>,
) -> APIResult<Json<IdentificationOutputData>> {
    let masses: Vec<String> = to_vec_string(data.masses)
        .map_err(|e| MetIDParseError::nums("mz", e))?;

    // An omitted `ccss` field becomes an empty vector.
    let mut ccss: Vec<String> = match data.ccss {
        Some(input) => to_vec_string(input)
            .map_err(|e| MetIDParseError::nums("ccs", e))?,
        None => Vec::new(),
    };

    let mass_error: String = to_num(data.mass_error)
        .map_err(|e| MetIDParseError::nums("mass_error", e))?
        .to_string();

    let mzwindow: f64 = to_num(data.mzwindow)
        .map_err(|e| MetIDParseError::nums("mz_window", e))?;

    let typ: Vec<String> = match data.typ {
        StringInput::OneStr(value) if value == "All" => vec![
            "Exogenous".to_string(),
            "Endogenous".to_string(),
            "Unspecified".to_string(),
        ],
        other => strnum_to_vec_string(other)
            .map_err(|e| MetIDParseError::string("typ", e))?,
    };

    let adducts: Vec<String> = if data.adducts == StringInput::OneStr("All".to_string()) {
        vec!["Phenolic Hydroxyls".to_string()]
    } else {
        strnum_to_vec_string(data.adducts)
            .map_err(|e| MetIDParseError::string("adducts", e))?
    };

    let requested_ccs =
        !ccss.is_empty() && ccss.first().is_some_and(|value| value != "-1");

    let ccs_included = requested_ccs && ccss.len() == masses.len();

    if requested_ccs && !ccs_included {
        println!(
            "The number of m/z values ({}) does not match the number of CCS values ({}); CCS will not be used for identification.",
            masses.len(),
            ccss.len()
        );
    }

    if !ccs_included {
        ccss = vec!["0".to_string(); masses.len()];
    }

    let processed = sql_handler(
        data.met,
        data.mat,
        typ,
        adducts,
        mass_error,
        masses,
        ccss,
        mzwindow,
        ccs_included,
    );

    Ok(Json(IdentificationOutputData(
        concatenate_maps_grouped(processed),
    )))
}

fn concatenate_maps_grouped(
    nested: Vec<Vec<HashMap<&'static str, String>>>,
) -> HashMap<&'static str, Vec<String>> {
    let mut result: HashMap<&'static str, Vec<String>> = HashMap::new();

    for inner_vec in nested {
        for map in inner_vec {
            for (key, value) in map {
                result.entry(key).or_default().push(value);
            }
        }
    }

    result
}

// Handler
pub async fn process_items(
    Json(data): Json<InputData>,
) -> Result<Json<OutputData>, JsonParseError> {
    let processed: Vec<f32> = data.0.into_iter().map(|value| value + 5.0).collect();
    Ok(Json(OutputData(processed)))
}