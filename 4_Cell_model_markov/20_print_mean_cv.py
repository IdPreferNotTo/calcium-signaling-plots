import numpy as np
import os
import warnings


def print_mean_CV_to_file():
    home = os.path.expanduser("~")

    N = 0
    n_cha = 5
    n_clu = 10
    taus = np.logspace(-1, 2, 50)
    js = np.logspace(-3, 0, 50)
    with open(home + "/Data/Calcium/data/markov_ca_mean_CV_N{:d}_n{:d}_no_adap.dat".format(n_clu, n_cha),
              "w") as outfile:
        outfile.write("# tau | j | <T> | CV | n \n")
        for i in range(50):
            for j in range(50):
                print(i, j)
                N += 1
                tau = taus[j]
                j = js[i]


                #file_str = home + "/CLionProjects/PhD/calcium_spikes_markov/out/spike_times_markov_taua{:.2e}_ampa{:.2e}_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(taua, ampa, tau, j, n_clu)
                file_str = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_no_adap/spike_times_markov_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu)

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
