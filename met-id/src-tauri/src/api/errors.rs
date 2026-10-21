use std::fmt;
use axum::{
    extract::rejection::JsonRejection, 
    http::StatusCode, 
    response::{IntoResponse, Response}, 
    Json
};
use serde_json::json;
use serde::Serialize;

use crate::api::parsing::{ParseNumsError, ParseStringError};



#[derive(Debug, thiserror::Error)]
pub enum MetIDError {

    #[error(transparent)]
    Json(#[from] JsonRejection),

    #[error(transparent)]
    Parsing(#[from] MetIDParseError),

    #[error(transparent)]
    Sql(#[from] rusqlite::Error),
}

#[derive(Serialize)]
struct ErrorBody {
    error: String
}

impl IntoResponse for MetIDError {
    fn into_response(self) -> Response {
        let status = match &self {
            MetIDError::Json(_) | MetIDError::Parsing(_)  => StatusCode::BAD_REQUEST,
            MetIDError::Sql(_) => StatusCode::INTERNAL_SERVER_ERROR,
        };
        let body = Json(ErrorBody {
            error: self.to_string()
        });
        (status, body).into_response()
    }
}

#[derive(Debug, thiserror::Error)]
pub enum MetIDParseError {
    #[error("invalid argument '{field}': {source}")]
    FieldNums {
        field: &'static str,
        #[source]
        source: ParseNumsError
    },

    #[error("invalid argument '{field}': {source}")]
    FieldString {
        field: &'static str,
        #[source]
        source: ParseStringError,
    },

    #[error("{0}")]
    Message(String)
}

impl MetIDParseError {
    pub fn nums(field: &'static str, source: ParseNumsError) -> Self {
        Self::FieldNums { field, source }
    }

    pub fn string(field: &'static str, source: ParseStringError) -> Self {
        Self::FieldString { field, source }
    }

    
}

pub type APIResult<T> = Result<T, MetIDError>;

pub struct JsonParseError(pub JsonRejection);

impl fmt::Debug for JsonParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.0)
    }
}

impl IntoResponse for JsonParseError {
    fn into_response(self) -> Response {
        match &self.0 {
            JsonRejection::MissingJsonContentType(_) => {
                (
                    StatusCode::BAD_REQUEST,
                    Json(json!({ "error": "Missing content-type: application/json" })),
                )
                    .into_response()
            }

            JsonRejection::JsonDataError(_) => {
                (
                    StatusCode::BAD_REQUEST,
                    Json(json!({
                        "error": "Invalid JSON Data"
                    })),
                )
                    .into_response()
            }

            _ => (
                StatusCode::UNPROCESSABLE_ENTITY,
                Json(json!({ "error": "Invalid request" })),
            )
                .into_response(),
        }
    }
}

impl From<JsonRejection> for JsonParseError {
    fn from(rejection: JsonRejection) -> Self {
        Self(rejection)
    }
}