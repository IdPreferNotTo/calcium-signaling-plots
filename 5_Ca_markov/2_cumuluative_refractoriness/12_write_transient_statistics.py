import os
import numpy as np
import warnings
from scipy.optimize import curve_fit

import functions as fc
import default_parameters as df

def print_mean_CV_with_adap_to_file():
    K = 10
    N = 5
    tau = 1
    p = 0.06
    tau_ers = np.logspace(1, 3, 41)
    eps_ers = np.logspace(-2, 0, 41)

    home = os.path.expanduser("~")
    with open(home + f"/Data/calcium_spikes_theory/markov_ca_transient_statistics_tau{tau:.2e}_j{p:.2e}_K{K:d}_N{N:d}_adap.dat",
              "w") as outfile:
        outfile.write("# tau_a | eps_a | NTr | T0 | T8 | var NTr | var T0 | var T8 | n \n")
        I = 0
        for tau_er in tau_ers:
            for eps_er in eps_ers:
                data_isi = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er)
                rows, cols = data_isi.shape
                idx_max = cols
                idxs = np.arange(idx_max)
                means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs]
                popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=(100, 150, 2))
                T0 = popt[0] - T0
                T8 = popt[1]
                NTr = popt[2]
                varT0 = pcov[0, 0]
                varT8 = pcov[1, 1]
                varNTr = pcov[1, 1]
                I += 1
                print(I)
            outfile.write(f"{tau_er:.2e} {eps_er:.2e} {NTr:.2e} {T0:.2e} {T8:.2e} {varT0:.2e} {varT8:.2e} {varNTr:.2e} \n")