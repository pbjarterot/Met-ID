use serde::{Deserialize, Serialize};
use super::parsing::NumInput;
use super::errors::APIResult;
use axum::Json;
use crate::api::parsing::to_num;
use crate::files::read_ms1_imzml;

#[derive(Debug, Deserialize)]
pub struct ImzmlInputData {
    path: String,
    n: NumInput,
    tolerance: NumInput
}

#[derive(Debug, Serialize)]
pub struct ImzmlOutputData {
    mz: Vec<f64>,
    intensity: Vec<f64>,
}

pub async fn read_imzml(Json(data): Json<ImzmlInputData>) 
    -> APIResult<Json<ImzmlOutputData>> {
    


    let res: Vec<Vec<f64>> = read_ms1_imzml(
        data.path, 
        to_num(data.n).unwrap() as usize, 
        to_num(data.tolerance).unwrap()
    );
    let (mz, int): (Vec<f64>, Vec<f64>) = res
    .into_iter()
    .map(|pair| (pair[0], pair[1]))
    .unzip();

    Ok(
        Json(ImzmlOutputData { 
            mz: mz, 
            intensity: int 
        })
    )
} 