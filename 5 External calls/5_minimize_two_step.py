import numpy as np
import subprocess
import os

from scipy.optimize import curve_fit
from scipy.optimize import minimize

import functions as fc


def get_non_adaptive_isis(xs):
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/langevin_sequence_fit_no_adap/" \
                       "cmake-build-release/langevin_sequence_fit_no_adap"
    result = subprocess.run([exetuable, f"{xs[0]:.8f}", f"{xs[1]:.8f}"],
                            stdout=subprocess.PIPE)
    result_string = result.stdout.decode('utf-8').strip()
    result_array = np.fromstring(result_string, dtype=float, sep=' ')
    return result_array


def func_non_adap(xs, args):
    t0_aim = args[0]
    cv_aim = args[1]
    isis = get_non_adaptive_isis(xs)
    t0 = np.mean(isis)
    cv = np.std(isis) / t0
    err = np.abs((t0 - t0_aim) / t0_aim) + np.abs((cv - cv_aim) / cv_aim)
    print(xs, "->", f"{t0:.4f}", f"{cv:.4f}", "->", err)
    return err


def get_adaptive_isis(xs, args):
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/langevin_transient_fit_with_adap/" \
                       "cmake-build-release/langevin_transient_fit_with_adap"
    result = subprocess.run([exetuable, f"{args[0]:.8f}", f"{args[1]:.8f}", f"{xs[0]:.8f}", f"{xs[1]:.8f}"],
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

    t0 = mean_ti[0]
    t8 = mean_ti[-1]
    func = lambda x, ntr: fc.exponential_Ti(x, t0, t8, ntr)  # fix t0, t8
    popt, pcov = curve_fit(func, idxs, mean_ti, p0=[2.])
    ntr = popt[0]
    return ntr, t0, t8, cv


def func_adap(xs, args):
    params = args[0]
    aim = args[1]

    isis = get_adaptive_isis(xs, params)
    ys = calculate_statistics(isis)
    ntr = ys[0]
    t0 = ys[1]
    t8 = ys[2]
    cv = ys[3]
    err = np.abs((ntr - aim[0]) / aim[0]) + np.abs((t8 - aim[1]) / aim[1])
    print(xs, "->", f"{ntr:.3f}", f"{t0:.2f}", f"{t8:.2f}", f"{cv:.4f}", " ->", err)
    return err


if __name__ == "__main__":
    aims_adap = np.asarray(
        [[1.1, 343], [4.6, 121], [7.4, 224], [2.1, 219], [3.5, 187], [2.4, 284],
         [0.9, 187], [2.4, 311], [1.6, 297], [6.4, 146], [3.3, 139], [3.8, 136],
         [4.9, 164], [1.8, 123], [1.6, 91], [1.8, 165], [3.8, 403], [2.6, 143],
         [4.3, 194], [4.2, 194], [1.9, 174], [7.4, 232], [4.7, 109], [5.8, 175]])

    params = np.asarray(
        [[10.2, 8.82e-03], [14.1, 9.90e-03], [11.6, 1.16e-02], [5.34, 2.14e-02], [4.20, 1.86e-02], [12.9, 1.19e-02],
         [2.71, 2.53e-02], [5.48, 1.47e-02], [5.48, 1.47e-02], [5.48, 1.47e-02], [52.1, 6.80e-03], [6.46, 1.51e-02],
         [11.4, 1.13e-02], [2.43, 3.13e-02], [12.2, 8.48e-03], [2.78, 2.67e-02], [3.93, 2.05e-02], [10.3, 8.65e-03],
         [18.6, 6.02e-03], [18.6, 6.02e-03], [7.56, 1.17e-02], [15.5, 7.47e-03], [7.71, 1.09e-02], [6.94, 1.74e-02]])

    x0s = np.asarray(
        [[4.14e+02, 1.78e-01], [6.49e+02, 5.5e-02], [1.94e+03, 4.53e-02], [5.12e+02, 1.44e-01], [6.50e+02, 3.68e-02], [8.02e+02, 1.60e-01],
         [1.75e+02, 8.52e-02], [7.59e+02, 6.30e-02], [5.01e+02, 9.45e-02], [9.81e+02, 1.90e-02], [6.00e+02, 1.15e-01], [5.79e+02, 5.23e-02],
         [9.36e+02, 5.84e-02], [2.48e+02, 6.63e-02], [1.54e+02, 8.71e-02], [3.10e+02, 6.25e-02], [1.60e+03, 4.10e-02], [4.15e+02, 5.27e-02],
         [9.31e+02, 4.84e-02], [9.40e+02, 4.83e-02], [3.65e+02, 8.35e-02], [1.92e+03, 3.41e-02], [5.37e+02, 2.38e-02], [1.19e+03, 5.05e-02]])


    a = (0.5 - 0.3) / (10 * 0.4623)
    x0 = np.asarray([10, 0.01]) #tau, p
    y_aim = [100, 0.20] # t0, cv

    res = minimize(func_non_adap, x0, args=y_aim, method="Nelder-Mead", options={'fatol': 0.05})
    # i = 1
    # print(params[i])
    # res = minimize(func_adap, x0s[i], args=[params[i], aims_adap[i]], method="Nelder-Mead",
    #                options={'fatol': 0.01, 'return_all': True})
    print(res)
