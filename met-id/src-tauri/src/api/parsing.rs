//use std::str::FromStr;
use serde::{Serialize, Deserialize};
use std::fmt;

use crate::api::errors::MetIDParseError;

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(untagged)]
pub enum NumInput {
    OneNum(f64),
    OneStr(String),
    ManyStr(Vec<String>),
    ManyNum(Vec<f64>),
}


#[derive(Debug, Clone, Deserialize, Serialize, PartialEq)]
#[serde(untagged)]
pub enum StringInput {
    OneStr(String),
    ManyStr(Vec<String>),
    OneNum(f64),
    ManyNum(Vec<f64>)
}

#[derive(Debug, thiserror::Error)]
pub enum ParseStringError {
    NotString(String),
    EmptyString(String)
}

impl fmt::Display for ParseStringError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ParseStringError::NotString(s) => write!(f, "Not a String: {}", s),
            ParseStringError::EmptyString(s) => write!(f, "String is empty: {s}"),
        }
    }
}

pub fn strnum_to_vec_string(input: impl Into<StringInput>) -> Result<Vec<String>, ParseStringError> {
    match input.into() {
        StringInput::OneNum(z) => Ok(vec![z.to_string()]),
        StringInput::ManyNum(z) => Ok(z.into_iter().map(|n| n.to_string()).collect()),
        StringInput::OneStr(a) => {
            match require_nonempty_str("", &a) {
                Ok(()) => Ok(vec![a]),
                Err(_e) => Err(ParseStringError::EmptyString("Empty string".to_string())),
            }
        },
        StringInput::ManyStr(z) => Ok(z)
    }
}

#[derive(Debug, thiserror::Error)]
pub enum ParseNumsError {
    ParseFloat(std::num::ParseFloatError),
    NotNumber(String),
    Empty
}

impl fmt::Display for ParseNumsError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ParseNumsError::ParseFloat(e) => write!(f, "{}", e),
            ParseNumsError::NotNumber(s) => write!(f, "Not a number: {}", s),
            ParseNumsError::Empty => write!(f, "Empty"),
        }
        
    }
}

impl fmt::Display for NumInput {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl From<f64> for NumInput {
    fn from(v: f64) -> Self { Self::OneNum(v) }
}
impl From<&str> for NumInput {
    fn from(v: &str) -> Self { Self::OneStr(v.to_string()) }
}
impl From<String> for NumInput {
    fn from(v: String) -> Self { Self::OneStr(v) }
}
impl From<Vec<String>> for NumInput {
    fn from(v: Vec<String>) -> Self { Self::ManyStr(v) }
}
impl From<Vec<f64>> for NumInput {
    fn from(v: Vec<f64>) -> Self { Self::ManyNum(v) }
}
impl<const N: usize> From<[f64; N]> for NumInput {
    fn from(v: [f64; N]) -> Self { Self::ManyNum(v.to_vec()) }
}

pub fn to_num(input: NumInput) -> Result<f64, ParseNumsError> {
    match input.into() {
        NumInput::OneNum(x) => Ok(x),
        NumInput::OneStr(x) => {
            match x.parse::<f64>() {
                Ok(a) => return Ok(a),
                Err(_) => return Err(ParseNumsError::NotNumber("Input cannot be parsed, input needs to be a single number".to_owned().to_string())),
            }
        },
        _ => Err(ParseNumsError::NotNumber("Input cannot be parsed, input needs to be a single number".to_string())),

    }
}

pub fn to_vec_string(input: NumInput) -> Result<Vec<String>, ParseNumsError> {
    match input.into() {
        NumInput::OneNum(x) => {
            Ok(vec![x.to_string()]) 
        },
        NumInput::ManyNum(xs) => {
            match require_nonempty_vec("", &xs) {
                Ok(()) => Ok(xs.into_iter().map(|v| v.to_string()).collect()),
                Err(_) => Err(ParseNumsError::Empty),
            }
            
        },
        NumInput::OneStr(s) => {
            match require_nonempty_str("", &s) {
                Ok(()) => Ok(vec![s]),
                Err(_) => Err(ParseNumsError::Empty),
            }
        },
        NumInput::ManyStr(ss) => {
            match require_nonempty_vec("", &ss) {
                Ok(()) => Ok(ss),
                Err(_) => Err(ParseNumsError::Empty),
            }
        }
        
    }
}


pub fn require_nonempty_str(field: &'static str, s: &str) -> Result<(), ParseStringError> {
    if s.trim().is_empty() {
        Err(ParseStringError::EmptyString(format!("'{field}' must not be empty")))
    } else {
        Ok(())
    }
}

pub fn require_nonempty_vec<T>(field: &'static str, v: &[T]) -> Result<(), MetIDParseError> {
    if v.is_empty() {
        Err(MetIDParseError::Message(format!("'{field}' must not be empty")))
    } else {
        Ok(())
    }
}
/* 
pub fn require_no_empty_items(field: &'static str, v: &[String]) -> Result<(), MetIDParseError> {
    if v.iter().any(|s| s.trim().is_empty()) {
        Err(MetIDParseError::Message(format!("'{field}' contains empty values")))
    } else {
        Ok(())
    }
}
*/