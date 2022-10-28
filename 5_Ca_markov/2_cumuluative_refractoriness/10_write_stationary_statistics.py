import os
import numpy as np
import warnings
from scipy.optimize import curve_fit

import functions as fc
import default_parameters as df

def print_mean_CV_with_adap_to_file():
    K = 10
    N = 5
    tau = 5
    p = 0.015
    tau_ers = np.logspace(1, 3, 41)
    eps_ers = np.logspace(-2, 0, 41)

    home = os.path.expanduser("~")
    with open(home + f"/Data/calcium_spikes_theory/markov_ca_stationary_statistics_tau{tau:.2e}_j{p:.2e}_K{K:d}_N{N:d}_adap.dat",
              "w") as outfile:
        outfile.write("# tau_a | eps_a | <T> | CV | n \n")
        I = 0
        for tau_er in tau_ers:
            for eps_er in eps_ers:
                data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
                if data_isi.size < 50:
                    mean = np.nan
                    cv = np.nan
                    n = 0
                else:
                    mean = np.mean(data_isi[50:])
                    std = np.std(data_isi[50:])
                    cv = std/mean
                    n = data_isi.size
                    I += 1
                print(I)
                outfile.write(f"{tau_er:.2e} {eps_er:.2e} {mean:.2e} {cv:.2e} {n:d} \n")

print_mean_CV_with_adap_to_file()