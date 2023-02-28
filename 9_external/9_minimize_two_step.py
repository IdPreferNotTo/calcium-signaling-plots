import numpy as np
import subprocess
import os

from scipy.optimize import curve_fit
from scipy.optimize import minimize

import functions as fc

def get_non_adaptive_isis(xs):
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin_fit/" \
               "cmake-build-release/calcium_spikes_langevin_fit"
    result = subprocess.run([exetuable, f"{xs[0]:.8f}", f"{xs[1]:.8f}"],
                    stdout=subprocess.PIPE)
    result_string = result.stdout.decode('utf-8').strip()
    result_array = np.fromstring(result_string, dtype=float, sep=' ')
    return result_array


def func_non_adap(xs):
    t0_aim = 16
    cv_aim = 0.27
    isis = get_non_adaptive_isis(xs)
    t0 = np.mean(isis)
    cv = np.std(isis)/t0
    err = np.abs((t0 - t0_aim)/t0_aim) + np.abs((cv- cv_aim)/cv_aim)
    print(xs, "->", f"{t0:.4f}", f"{cv:.4f}", "->", err)
    return err


def get_adaptive_isis(xs, args):
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin_transient/" \
               "cmake-build-release/calcium_spikes_langevin_transient"
    result = subprocess.run([exetuable, f"{args[0]:.8f}", f"{args[1]:.8f}", f"{xs[0]:.8f}", f"{xs[1]:.8f}"],
                    stdout=subprocess.PIPE)
    result_string = result.stdout.decode('utf-8').strip().split("\n")
    result_string = [str.split() for str in result_string]
    result_array = np.array(result_string, float)
    return result_array


def calculate_statistics(isis):
    # Calculate CV from stationary interspike intervals (Ti with i 20 - 25)
    t0 = np.mean(isis[:, 0])
    t8 = np.mean(isis[:, -1])

    stationary_ti = isis[:, -5:-1]
    stationary_ti = stationary_ti.flatten()
    cv = np.std(stationary_ti)/np.mean(stationary_ti)

    # Calculate ntr, t0, t8 from <Ti> with exponential fit
    mean_ti= np.mean(isis, axis=0)
    idxs = np.arange(mean_ti.size)

    func = lambda x, ntr: fc.exponential_Ti(x, t0, t8, ntr) # fix t0, t8
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
    err = np.abs((ntr - aim[0])/aim[0]) + np.abs((t8 - aim[1])/aim[1])
    print(xs, "->", f"{ntr:.3f}", f"{t0:.2f}", f"{t8:.2f}", f"{cv:.4f}", " ->", err)
    return err



if __name__ == "__main__":
    params = np.asarray([[10.28, 8.779e-03], [13.38, 9.923e-03], [11.37, 1.168e-02],
                       [4.532, 2.288e-02], [4.180, 1.866e-02], [9.692, 1.308e-02],
                       [2.741, 2.510e-02], [5.492, 1.471e-02], [5.492, 1.471e-02],
                       [5.492, 1.471e-02], [65, 6.7e-03], [6.537, 1.501e-02],
                       [11.27, 1.139e-02], [2.407, 3.133e-02], [12.08, 8.519e-03],
                       [2.781, 2.666e-02], [3.988, 2.029e-02], [10.27, 8.700e-03],
                       [18.74, 6.000e-03], [18.74, 6.000e-03], [7.619, 1.163e-02],
                       [15.49, 7.463e-03], [7.777, 1.087e-02], [6.850, 1.749e-02]])
    aims = np.asarray([[1.1, 343], [4.6, 121], [7.4, 224],
                       [2.1, 219], [3.5, 187], [2.4, 284],
                       [0.9, 187], [2.4, 311], [1.6, 297],
                       [6.4, 146], [3.3, 139], [3.8, 136],
                       [4.9, 164], [1.8, 123], [1.6, 91],
                       [1.8, 165], [3.8, 403], [2.6, 143],
                       [4.3, 194], [4.2, 194], [1.9, 174],
                       [7.4, 232], [4.7, 109], [5.8, 175]])

    x0 = np.asarray([2.4, 0.031]) # tau_i, p, tau_er, eps

    res = minimize(func_non_adap, x0, method="Nelder-Mead", options={'fatol': 0.01})
    # i = 10
    # print(params[i])
    # res = minimize(func_adap, x0, args=[params[i], aims[i]], method="Nelder-Mead", options={'fatol': 0.01, 'return_all': True})
    print(res)