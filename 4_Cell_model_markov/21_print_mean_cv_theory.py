import numpy as np
import os
import warnings

import functions as fc

def print_mean_CV_theory_to_file():
    home = os.path.expanduser("~")

    N = 0
    n_cha = 5
    n_clu = 10
    taus = np.logspace(-1, 2, 50)
    js = np.logspace(-3, 0, 50)
    with open(home + "/Data/Calcium/data/langevin_theory_ca_mean_CV_N{:d}_n{:d}_no_adap.dat".format(n_clu, n_cha),
              "w") as outfile:
        outfile.write("# tau | j | <T> | CV | n \n")
        for i in range(50):
            for j in range(50):
                print(i, j)
                N += 1
                tau = taus[j]
                j = js[i]
                r0 = 1


                outfile.write("{:.2e} {:.2e} {:.2e} {:.2e} \n".format(tau, j, T, CV))


def find_j_tau_for(ISI, CV):
    home = os.path.expanduser("~")
    data = np.loadtxt(home + "/Data/Calcium/data/langevin_ca_mean_CV_N10_n5_no_adap.dat")
    print(data)

if __name__ == "__main__":
    find_j_tau_for(1, 1)
