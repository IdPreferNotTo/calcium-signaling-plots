import numpy as np
import os
from scipy.optimize import curve_fit

import styles as st
import functions as fc

def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))


if __name__ == "__main__":
    home = os.path.expanduser("~")

    N = 0
    n_cha = 5
    n_clu = 10
    taua = 500
    ampa = 0.05
    folder = "/Data/Calcium/data"
    with open(home + folder + f"/markov_ca_mean_CV_N{N:d}_n{n_cha:d}_fix_adap_taua{taua:.2e}_ampa{ampa:.2e}.dat", "w") as outfile:
        for j in np.logspace(-3, 0, 31):
            for tau in np.logspace(-1, 2, 31):
                N += 1
                file_str = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap_fix_adap_para/spike_times_markov_ip1.00_taua5.00e+02_ampa5.00e-02_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(
                    tau, j, n_clu)