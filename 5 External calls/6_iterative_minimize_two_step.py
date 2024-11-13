import numpy as np
import subprocess
import os

from scipy.optimize import curve_fit
from scipy.optimize import minimize

import functions as fc


def get_renewal_isis(xs):
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/langevin_sequence_fit_no_adap/" \
                       "cmake-build-release/langevin_sequence_fit_no_adap"
    result = subprocess.run([exetuable, f"{xs[0]:.8f}", f"{xs[1]:.8f}"],
                            stdout=subprocess.PIPE)
    result_string = result.stdout.decode('utf-8').strip()
    result_array = np.fromstring(result_string, dtype=float, sep=' ')
    return result_array


def loss_function_renewal(paras, stats):
    t0_aim = stats[0]
    cv_aim = stats[1]
    isis = get_renewal_isis(paras)
    t0 = np.mean(isis)
    cv = np.std(isis) / t0
    err = np.abs((t0 - t0_aim) / t0_aim) + np.abs((cv - cv_aim) / cv_aim)
    print(paras, "->", f"{t0:.4f}", f"{cv:.4f}", "->", err)
    return err


def get_nonrenewal_isis(pars_ci, pars_cer):
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/langevin_transient_fit_with_adap/" \
                       "cmake-build-release/langevin_transient_fit_with_adap"
    result = subprocess.run([exetuable, f"{pars_ci[0]:.8f}", f"{pars_ci[1]:.8f}", f"{pars_cer[0]:.8f}", f"{pars_cer[1]:.8f}"],
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
    return t0, cv, ntr, t8


def loss_function_nonrenewal(xs, args):
    pars_ci = args[0]
    stats_cer = args[1]
    isis = get_nonrenewal_isis(pars_ci, xs)
    ys = calculate_statistics(isis)
    T0 = ys[0]
    CV = ys[1]
    Ntr = ys[2]
    T8 = ys[3]
    err =  np.abs(Ntr - stats_cer[0])/stats_cer[0]  + np.abs(T8 - stats_cer[1])/stats_cer[1]
    print(xs, "->", f"{T0:.2f}", f"{CV:.2f}", f"{Ntr:.2f}", f"{T8:.2f}", " ->", err)
    return err


if __name__=="__main__":
    stats_all = np.asarray(
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
         [41, 0.12, 2.6, 143],
         [46, 0.09, 4.3, 194],
         [46, 0.09, 4.2, 194],
         [31, 0.14, 1.9, 174],
         [36, 0.10, 7.4, 232],
         [36, 0.14, 4.7, 109],
         [15, 0.15, 5.8, 175]])

    pars_all_it0 = np.asarray(
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

    pars_all_it1 = np.asarray(
        [[6.09, 1.24e-02, 4.20e+02, 1.09e-01],
         [61.7, 7.30e-03, 7.36e+02, 8.73e-02],
         [99.3, 7.87e-03, 2.46e+03, 8.32e-02],
         [3.22, 2.73e-02, 5.18e+02, 9.78e-02],
         [7.11, 1.32e-02, 7.18e+02, 5.56e-02],
         [9.65, 1.31e-02, 8.00e+02, 1.31e-01],
         [2.27, 2.92e-02, 1.75e+02, 6.94e-02],
         [5.95, 1.40e-02, 7.87e+02, 6.71e-02],
         [4.42, 1.72e-02, 4.83e+02, 8.02e-02],
         [15.6, 8.47e-03, 1.06e+03, 4.14e-02],
         [53.3, 6.89e-03, 6.04e+02, 1.17e-01],
         [18.3, 1.00e-02, 6.48e+02, 9.02e-02],
         [115., 7.52e-03, 1.17e+03, 1.05e-01],
         [2.21, 3.34e-02, 2.29e+02, 5.99e-02],
         [10.2, 9.23e-03, 1.56e+02, 7.55e-02],
         [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
         [7.24, 1.42e-02, 1.69e+03, 6.79e-02],
         [13.9, 7.35e-03, 3.82e+02, 6.55e-02],
         [28.5, 5.03e-03, 9.54e+02, 5.86e-02],
         [35.4, 4.72e-03, 9.65e+02, 6.54e-02],
         [8.56, 1.08e-02, 3.53e+02, 9.42e-02],
         [77.1, 4.93e-03, 2.32e+03, 6.04e-02],
         [12.4, 8.30e-03, 5.46e+02, 3.36e-02],
         [73.3, 1.11e-02, 1.53e+03, 1.08e-01]])

    pars_all_it2 = np.asarray(
        [[6.09, 1.24e-02, 4.20e+02, 1.09e-01],
         [19.5, 8.95e-03, 7.06e+02, 6.44e-02],
         [22.4, 9.54e-03, 2.10e+03, 6.01e-02],
         [3.22, 2.73e-02, 5.18e+02, 9.78e-02],
         [7.11, 1.32e-02, 7.18e+02, 5.56e-02],
         [9.65, 1.31e-02, 8.00e+02, 1.31e-01],
         [1.84, 3.53e-02, 1.83e+02, 6.42e-02],
         [5.95, 1.40e-02, 7.87e+02, 6.71e-02],
         [4.42, 1.72e-02, 4.83e+02, 8.02e-02],
         [9.92, 1.04e-02, 1.04e+03, 3.13e-02],
         [53.3, 6.89e-03, 6.04e+02, 1.17e-01],
         [9.39, 1.26e-02, 5.91e+02, 6.72e-02],
         [23.5, 9.06e-03, 1.03e+03, 8.33e-02],
         [2.41, 3.14e-02, 2.28e+02, 6.75e-02],
         [10.2, 9.23e-03, 1.56e+02, 7.55e-02],
         [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
         [7.24, 1.42e-02, 1.69e+03, 6.79e-02],
         [13.9, 7.35e-03, 3.82e+02, 6.55e-02],
         [28.5, 5.03e-03, 9.54e+02, 5.86e-02],
         [35.4, 4.72e-03, 9.65e+02, 6.54e-02],
         [8.56, 1.08e-02, 3.53e+02, 9.42e-02],
         [32.3, 5.80e-03, 2.09e+03, 4.88e-02],
         [12.4, 8.30e-03, 5.46e+02, 3.36e-02],
         [77.0, 1.11e-02, 1.55e+03, 1.09e-01]])

    cvs_hek = stats_all[:, 1]
    cvs_model_it0 = [0.07, 0.14, 0.23, 0.14, 0.28, 0.11, 0.22, 0.18, 0.14, 0.28, 0.08, 0.21, 0.18, 0.25, 0.10, 0.26, 0.28, 0.14, 0.11, 0.11, 0.15, 0.16, 0.18, 0.26]
    cvs_model_it1 = [0.11, 0.09, 0.10, 0.18, 0.20, 0.13, 0.26, 0.18, 0.16, 0.14, 0.09, 0.12, 0.09, 0.29, 0.11, 0.26, 0.21, 0.11, 0.09, 0.08, 0.14, 0.08, 0.13, 0.11]
    cvs_aim_it0 = stats_all[:, 1]
    cvs_aim_it1 = cvs_aim_it0 * (1. + (cvs_hek - cvs_model_it0)/np.maximum(cvs_hek, cvs_model_it0))
    cvs_aim_it2 = cvs_aim_it1 * (1. + (cvs_hek - cvs_model_it1)/np.maximum(cvs_hek, cvs_model_it1))
    stats_all[:, 1] = cvs_aim_it2
    idx_of_interest = [1, 2, 6, 9, 11, 12, 13, 21, 23]
    idx = 2
    stats = stats_all[idx]
    pars = pars_all_it1[idx]
    stats_ci = stats[:2]
    stats_cer = stats[2:]
    pars_ci = [22.4, 9.54e-03] #pars[:2]
    pars_cer = [2.36e+03, 5.39e-02] #pars[2:]

    if False:
        res = minimize(loss_function_renewal, pars_ci, args=stats_ci, method="Nelder-Mead",
                   options={'xatol':0.02,'fatol': 0.02, 'disp':True})

        pars_ci = res["x"]

    if True:
        res = minimize(loss_function_nonrenewal, pars_cer, args=[pars_ci, stats_cer], method="Nelder-Mead",
                       options={'fatol': 0.05, 'adaptive': True})

        print(res["x"])
        pars_cer = res["x"]

    isis = get_nonrenewal_isis(pars_ci, pars_cer)
    out_stats = calculate_statistics(isis)
    print(f"{out_stats[0]:.0f}, {out_stats[2]:.1f}, {out_stats[3]:.0f}, {out_stats[1]:.2f}")