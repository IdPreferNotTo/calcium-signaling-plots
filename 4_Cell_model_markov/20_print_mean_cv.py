import numpy as np
import os
import warnings


def print_mean_CV_to_file():
    home = os.path.expanduser("~")

    K = 10
    N = 5
    with open(home + "/Data/calcium_spikes_theory/langevin_haenggi_ca_mean_CV_K{:d}_N{:d}_no_adap.dat".format(K, N),
              "w") as outfile:
        outfile.write("# tau | j | <T> | CV | n \n")
        for j in np.logspace(-3, 0, 61):
            for tau in np.logspace(-1, 2, 61):
                N += 1
                file_str = home + "/Data/calcium_spikes_langevin_haenggi/Data_no_adap/spike_times_langevin_ip1.00_tau{:.2e}_j{:.2e}_K{:d}_0.dat".format(tau, j, K)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ISIs = np.loadtxt(file_str, unpack=True)
                    if ISIs.size == 0:
                        T = np.infty
                        CV = 1.
                        n = ISIs.size
                    else:
                        T = np.mean(ISIs)
                        CV = np.std(ISIs)/np.mean(ISIs)
                        n = ISIs.size
                outfile.write("{:.2e} {:.2e} {:.2e} {:.2e} {:d}\n".format(tau, j, T, CV, n))


if __name__ == "__main__":
    print_mean_CV_to_file()
