extern crate nlopt;

use log::info;
use nlopt::{Algorithm, Nlopt};

fn lq(parameters: &[f64], x: f64) -> f64 {
    parameters[0] * x.powi(2) + parameters[1] * x + parameters[2]
}

fn objective(parameters: &[f64], x: &[f64], y: &[f64]) -> f64 {
    let mut err = vec![0.0; x.len()];
    for i in 0..x.len() {
        err[i] = y[i] - lq(parameters, x[i]);
    }

    let mut sum_err_sq: f64 = 0.0;
    for i in 0..x.len() {
        sum_err_sq += err[i].powi(2);
    }
    sum_err_sq
}

fn fit_crv(x0: Vec<f64>, data: &[(f64, f64)], algo: Algorithm) -> Vec<f64> {
    let x_array: Vec<f64> = data.iter().map(|&(x, _)| x).collect();
    let y_array: Vec<f64> = data.iter().map(|&(_, y)| y).collect();
    let objective_closure =
        |parameters: &[f64], _grad: Option<&mut [f64]>, _data: &mut &[(f64, f64)]| {
            objective(parameters, &x_array, &y_array) //, &[data[0].0], &[data[0].1])
        };
    let mut opt = Nlopt::new(algo, 3, objective_closure, nlopt::Target::Minimize, data);
    opt.set_xtol_abs(&[1e-60]).unwrap();
    opt.set_maxtime(5.0).unwrap();
    opt.set_maxeval(1_000_000).unwrap();

    let mut plsq: Vec<f64> = x0.clone();
    match opt.optimize(&mut plsq) {
        Ok(_) => {}
        Err(err) => {
            info!("Either the arguments are 0, 0, 0 or something is wrong, if the latter please contact patrik.bjarterot@uu.se");
            println!("something happened: {:?}", err);
            return vec![0.0, 0.0, 0.0];
        }
    }
    plsq
}

fn fit_quadratic_curve(data: &[(f64, f64)]) -> Vec<f64> {
    let mut plsq: Vec<f64> = vec![0.0, 0.0, 0.0];
    plsq = fit_crv(plsq, data, Algorithm::Newuoa);
    plsq = fit_crv(plsq, data, Algorithm::Bobyqa);
    plsq
}

#[tauri::command]
pub fn mass_error_regression(data: Vec<(f64, f64)>) -> (f64, f64, f64) {
    let fit_res: Vec<f64> = fit_quadratic_curve(&data);
    (fit_res[0], fit_res[1], fit_res[2])
}
