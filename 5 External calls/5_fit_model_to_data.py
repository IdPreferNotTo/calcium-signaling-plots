import numpy as np
import subprocess
import os

from scipy.optimize import curve_fit
from scipy.optimize import minimize

import functions as fc


def get_adaptive_isis(xs):
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/langevin_transient_fit_with_adap/" \
                       "cmake-build-release/langevin_transient_fit_with_adap"
    result = subprocess.run([exetuable, f"{xs[0]:.8f}", f"{xs[1]:.8f}", f"{xs[2]:.8f}", f"{xs[3]:.8f}"],
                            stdout=subprocess.PIPE)
    result_string = result.stdout.decode('utf-8').strip().split("\n")
    result_string = [str.split() for str in result_string]
    result_array = np.array(result_string, float)
    return result_array


def calculate_statistics(isis):
    # Calculate CV from stationary interspike intervals (Ti with i 20 - 25)
    stationary_ti = isis[:, -5:-1]
    stationary_ti = stationary_ti.flatten()
    cv = np.std(stationary_ti) / np.mean(stationary_ti)

    # Calculate ntr, t0, t8 from <Ti> with exponential fit
    mean_ti = np.mean(isis, axis=0)
    idxs = np.arange(mean_ti.size)
    popt, pcov = curve_fit(fc.exponential_Ti, idxs, mean_ti, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
    T0 = popt[0]
    T8 = popt[1]
    n_tr = popt[2]
    return T0, cv, n_tr, T8


def func_adap(xs, args):
    aim = args

    isis = get_adaptive_isis(xs)
    ys = calculate_statistics(isis)
    T0 = ys[0]
    CV = ys[1]
    Ntr = ys[2]
    T8 = ys[3]
    err =  np.abs(T0 - args[0])/args[0] + np.abs(CV - args[1])/args[1]  + np.abs(Ntr - args[2])/args[2]  + np.abs(T8 - args[3])/args[3]
    print(xs, "->", f"{T0:.2f}", f"{CV:.2f}", f"{Ntr:.2f}", f"{T8:.2f}", " ->", err)
    return err


if __name__ == "__main__":
    #t0, cv, ntr, t8
    aims = np.asarray(
        [[40, 0.12, 1.1, 343],
         [24, 0.11, 4.6, 121],
         [21, 0.12, 7.4, 224],
         [13, 0.18, 2.1, 219],
         [25, 0.20, 3.5, 187],
         [19, 0.12, 2.4, 284],
         [29, 0.31, 0.9, 187],
         [29, 0.17, 2.4, 311],
         [29, 0.17, 1.6, 297],
         [29, 0.17, 6.4, 146],
         [26, 0.07, 3.3, 139],
         [21, 0.15, 3.8, 136],
         [22, 0.12, 4.9, 164],
         [16, 0.27, 1.8, 123],
         [35, 0.11, 1.6, 91],
         [20, 0.26, 1.8, 165],
         [21, 0.20, 3.8, 403],
         [41, 0.1, 2.6, 143],
         [46, 0.09, 4.3, 194],
         [46, 0.09, 4.2, 194],
         [31, 0.14, 1.9, 174],
         [36, 0.10, 7.4, 232],
         [36, 0.14, 4.7, 109],
         [15, 0.15, 5.8, 175]])

    # tau, p, tau_er, eps_er
    x0s = np.asarray(
        [[10.2, 8.82e-03, 4.14e+02, 1.78e-01],
         [14.1, 9.90e-03, 6.49e+02, 5.5e-02],
         [11.6, 1.16e-02, 1.94e+03, 4.53e-02],
         [5.34, 2.14e-02, 5.12e+02, 1.44e-01],
         [4.20, 1.86e-02, 6.50e+02, 3.68e-02],
         [12.9, 1.19e-02, 8.02e+02, 1.60e-01],
         [2.71, 2.53e-02, 1.75e+02, 8.52e-02],
         [5.48, 1.47e-02, 7.59e+02, 6.30e-02],
         [5.48, 1.47e-02, 5.01e+02, 9.45e-02],
         [5.48, 1.47e-02, 9.81e+02, 1.90e-02],
         [52.1, 6.80e-03, 6.00e+02, 1.15e-01],
         [6.46, 1.51e-02, 5.79e+02, 5.23e-02],
         [11.4, 1.13e-02, 9.36e+02, 5.84e-02],
         [2.43, 3.13e-02, 2.48e+02, 6.63e-02],
         [12.2, 8.48e-03, 1.54e+02, 8.71e-02],
         [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
         [3.93, 2.05e-02, 1.60e+03, 4.10e-02],
         [10.3, 8.65e-03, 4.15e+02, 5.27e-02],
         [18.6, 6.02e-03, 9.31e+02, 4.84e-02],
         [18.6, 6.02e-03, 9.40e+02, 4.83e-02],
         [7.56, 1.17e-02, 3.65e+02, 8.35e-02],
         [15.5, 7.47e-03, 1.92e+03, 3.41e-02],
         [7.71, 1.09e-02, 5.37e+02, 2.38e-02],
         [6.94, 1.74e-02, 1.19e+03, 5.05e-02]])


    i = 3
    res = minimize(func_adap, x0s[i], args=aims[i], method="Nelder-Mead",
                   options={'fatol': 0.05, 'return_all': True})
    print(res)
