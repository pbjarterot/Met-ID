
use pyo3::{prelude::*, prepare_freethreaded_python};
use pyo3::types::{PyDict, PyString};

fn fit_quat(_data: &Vec<(f64, f64)>) -> PyResult<(f64, f64, f64)> {
    prepare_freethreaded_python();
    let mut val: (f64, f64, f64) = (0.0, 0.0, 0.0);
    Python::with_gil(|py| {

        let locals: &PyDict = PyDict::new(py);
        locals.set_item("code", PyString::new(py, r#"
from scipy.optimize import fmin
import numpy as np

def lq(parameters, x):
    return (parameters[0]*x**2) + parameters[1]*x + parameters[2]

def objective(parameters, x, y, w):
    err = y - lq(parameters, x)
    return np.sum(w*err**2)

def fit_quadratic_curve(data):
    x_array = np.array([x for x, y in data])
    y_array = np.array([y for x, y in data])
    x0 = [0, 0, 0]
    w = [1 for i in list(x_array)]
    plsq = fmin(objective, x0, args=(x_array, y_array, w), ftol=1e-5, maxiter=10000, disp=False)

    return plsq

result = tuple(fit_quadratic_curve(data))"#)).unwrap();
        locals.set_item("data", _data).unwrap();

        //execute the function
        if let Err(e) = py.run("exec(code)", Some(locals), None) {
            e.print(py);
        }

        // Retrieve the result from the locals dictionary
        let result: &PyAny = locals.get_item("result").unwrap();
        val = result.extract().unwrap();
    });
    Ok(val)

}


#[tauri::command]
pub fn mass_error_regression(data: Vec<(f64, f64)>) -> (f64, f64, f64) {
    let res: (f64, f64, f64) = fit_quat(&data).unwrap();

    (res.0, res.1, res.2) 
}