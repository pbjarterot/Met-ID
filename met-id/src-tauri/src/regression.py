import sys
from scipy.optimize import fmin
import numpy as np

def lq(parameters, x):
    return (parameters[0]*x**2) + parameters[1]*x + parameters[2]

def objective(parameters, x, y, w):
    err = y - lq(parameters, x)
    return np.sum(w*err**2)

def fit_quadratic_curve(data):
    x_array = np.array(data[0::2])
    y_array = np.array(data[1::2])
    x0 = [0, 0, 0]
    w = [1 for i in x_array]
    plsq = fmin(objective, x0, args=(x_array, y_array, w), ftol=1e-5, maxiter=10000, disp=False)

    return plsq

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) % 2 != 0:
        print("Usage: python script.py x1 y1 x2 y2 ... xn yn")
        sys.exit(1)

    try:
        data = [float(arg) for arg in sys.argv[1:]]
        result = tuple(fit_quadratic_curve(data))
        print(result)
    except ValueError:
        print("Invalid data format in the command line arguments.")