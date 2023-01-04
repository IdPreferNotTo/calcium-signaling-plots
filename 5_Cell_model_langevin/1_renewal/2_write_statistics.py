import numpy as np
import os
import warnings

import default_parameters as df

def print_mean_CV_to_file():
    home = os.path.expanduser("~")
    K = 10
    N = 5
    with open(home + "/Data/calcium_spikes_theory/langevin_isi_mean_CV_K{:d}_N{:d}_no_adap.dat".format(K, N),
              "w") as outfile:
        outfile.write("# tau | j | <T> | CV | n \n")
        print(np.logspace(0, 2, 41))
        for p in np.logspace(-3, -1, 41):
            for tau in np.logspace(0, 2, 41):
                N += 1
                print(N, p, tau)
                data_isi = df.load_spike_times_langevin(tau, p, cer=False)
                if data_isi.size < 2:
                    T = np.infty
                    CV = 1.
                    n = data_isi.size
                else:
                    T = np.mean(data_isi)
                    CV = np.std(data_isi)/np.mean(data_isi)
                    n = data_isi.size
                outfile.write("{:.2e} {:.2e} {:.2e} {:.2e} {:d}\n".format(tau, p, T, CV, n))


if __name__ == "__main__":
    print_mean_CV_to_file()