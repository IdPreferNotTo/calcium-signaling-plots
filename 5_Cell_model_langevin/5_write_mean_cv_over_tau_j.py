import numpy as np
import os
import warnings


def print_mean_CV_to_file():
    home = os.path.expanduser("~")

    variances = []
    rates = []
    num_spikes = []
    N = 0
    taua = 100
    ampa = 0.2
    N = 5
    K = 10
    taus = np.logspace(-1, 2, 61)
    js = np.logspace(-3, 0, 61)
    with open(home + "/Data/calcium_spikes_theory/langevin_strat_ca_mean_CV_K{:d}_N{:d}_no_adap.dat".format(K, N),
              "w") as outfile:
        outfile.write("# tau | j | 1/T | CV | n \n")
        for tau in taus:
            for j in js:
                print(j, tau)
                N += 1

                #ile_str = home + "/Data/calcium_spikes_langevin/spike_times_langevin_ip1.00_taua{:.2e}_ampa{:.2e}_tau{:.2e}_j{:.2e}_K{:d}_5.dat".format(taua, ampa, tau, j, K)
                file_str = home + "/Data/calcium_spikes_langevin_strat/Data_no_adap/spike_times_langevin_ip1.00_tau{:.2e}_j{:.2e}_K{:d}_0.dat".format(tau, j, K)

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ISIs = np.loadtxt(file_str, unpack=True)
                    if ISIs.size == 0:
                        r = 0
                        CV = 1.
                        n = ISIs.size

                    else:
                        r = 1./np.mean(ISIs)
                        CV = np.std(ISIs)/np.mean(ISIs)
                        n = ISIs.size

                outfile.write("{:.2e} {:.2e} {:.2e} {:.2e} {:d}\n".format(tau, j, r, CV, n))


if __name__ == "__main__":
    print_mean_CV_to_file()
