import numpy as np
import os
import warnings


def print_mean_CV_to_file():
    home = os.path.expanduser("~")

    N = 0
    n_cha = 5
    n_clu = 10
    taua = 500
    ampa = 0.05
    with open(home + "/Data/Calcium/data/markov_ca_mean_CV_N{:d}_n{:d}_fix_adap_taua{:.2e}_ampa{:.2e}.dat".format(n_clu, n_cha, taua, ampa),
              "w") as outfile:
        outfile.write("# tau | j | <T> | CV | n \n")
        for j in np.logspace(-3, 0, 31):
            for tau in np.logspace(-1, 2, 31):
                N += 1
                file_str = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap_fix_adap_para/spike_times_markov_ip1.00_taua5.00e+02_ampa5.00e-02_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu)
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
