use axum::{
    extract::State as AxumState,
    http::StatusCode,
    Json,
};
use serde::{Deserialize, Serialize};

use crate::sidecar::{
    predict_with_sidecar,
    CcsSidecarState,
};
use std::sync::Arc;

#[derive(Clone)]
pub struct ApiState {
    pub ccs: Arc<CcsSidecarState>,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum OneOrMany<T> {
    One(T),
    Many(Vec<T>),
}

impl<T> OneOrMany<T> {
    fn into_vec(self) -> Vec<T> {
        match self {
            Self::One(value) => vec![value],
            Self::Many(values) => values,
        }
    }
}

#[derive(Debug, Deserialize)]
pub(crate) struct CcsApiInput {
    smiles: OneOrMany<String>,

    // Accept either "adduct" or "adducts" in the JSON request.
    #[serde(default, alias = "adducts")]
    adduct: Option<OneOrMany<String>>,
}

#[derive(Debug, Serialize)]
pub(crate) struct CcsApiResult {
    smiles: String,
    adduct: String,
    predicted_ccs: f64,
    mz: f64,
}

pub async fn predict_ccs_api(
    AxumState(state): AxumState<ApiState>,
    Json(input): Json<CcsApiInput>,
) -> Result<Json<Vec<CcsApiResult>>, (StatusCode, String)> {
    let mut smiles = input.smiles.into_vec();

    let mut adducts = input
        .adduct
        .map(OneOrMany::into_vec)
        .unwrap_or_else(|| vec!["[M+H]+".to_string()]);

    if smiles.is_empty() || smiles.iter().any(|value| value.trim().is_empty()) {
        return Err((
            StatusCode::BAD_REQUEST,
            "Provide at least one non-empty SMILES".to_string(),
        ));
    }

    if adducts.is_empty() {
        adducts.push("[M+H]+".to_string());
    }

    // Use the default adduct for blank values.
    for adduct in &mut adducts {
        if adduct.trim().is_empty() {
            *adduct = "[M+H]+".to_string();
        }
    }

    // Allow one SMILES with multiple adducts, or multiple SMILES with one
    // adduct. Otherwise the lists must have matching lengths.
    if smiles.len() == 1 && adducts.len() > 1 {
        smiles.resize(adducts.len(), smiles[0].clone());
    } else if adducts.len() == 1 && smiles.len() > 1 {
        adducts.resize(smiles.len(), adducts[0].clone());
    }

    if smiles.len() != adducts.len() {
        return Err((
            StatusCode::BAD_REQUEST,
            "Provide one SMILES/adduct pair, or lists with matching lengths"
                .to_string(),
        ));
    }

    let mut results = Vec::with_capacity(smiles.len());

    for (smiles, adduct) in smiles.into_iter().zip(adducts) {
        let prediction = predict_with_sidecar(
            state.ccs.as_ref(),
            smiles.clone(),
            Some(adduct.clone()),
        )
        .await
        .map_err(|error| {
            let message = error.to_string();

            let status = if message.starts_with("Invalid SMILES:")
                || message.starts_with("Invalid adduct ")
            {
                StatusCode::BAD_REQUEST
            } else {
                StatusCode::BAD_GATEWAY
            };

            (status, message)
        })?;

        results.push(CcsApiResult {
            smiles,
            adduct,
            predicted_ccs: prediction.predicted_ccs,
            mz: prediction.mz,
        });
    }

    Ok(Json(results))
}