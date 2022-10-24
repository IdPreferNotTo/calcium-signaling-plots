import numpy as np
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    data = np.loadtxt(home + "/Data/calcium_spikes_theory/markov_ca_mean_CV_K10_N5_no_adap.dat")

    T_aim = 157
    CV_aim = 0.15
    d = 100
    for set in data:
        tau = set[0]
        j = set[1]
        T = set[2]
        CV = set[3]
        dT = abs(T_aim - T) / T_aim
        dCV = abs(CV_aim - CV) / CV_aim
        if (dT + dCV) / 2. < d:
            T0 = T
            CV0 = CV
            tau0 = tau
            j0 = j
            d = (dT + dCV) / 2.

    print(T0, CV0, tau0, j0)
