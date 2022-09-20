import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    home = os.path.expanduser("~")
    taus = [1.78, 5.62]
    js = [0.0355, 0.0126]
    #taus = [4.47, 1.26]
    #js = [0.0178, 0.0562]
    IP3s = [1.] #[1.2, 1.4, 1.8]
    for tau, j in zip(taus, js):
        print(tau, j)
        for IP3 in IP3s:
            folder_markov_ip3 = home + "/Data/calcium_spikes_markov/Data_no_adap_ip3/"
            p_isis = []
            file_markov = f"spike_times_markov_ip{IP3:.2f}_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
            warnings.simplefilter("ignore")
            isis_markov = np.loadtxt(folder_markov_ip3 + file_markov)

            with open(home + f"/Data/calcium_spikes_theory/Fano_factor_ip{IP3:.2f}_tau{tau:.2e}_j{j:.2e}.dat",
                      "w") as file:
                file.write("#Delta T mean-ISI std-ISI Fano-ISI \n")
                dts = np.logspace(-1, 3, 1000)
                mean_counts = []
                var_counts = []
                for dt in dts:
                    print(dt)
                    mean_count, var_count = fc.fano_factor_interspike_intervals(isis_markov, dt)
                    file.write(f"{dt:.2e} {mean_count:.2e} {var_count:.2e} {var_count/mean_count:.2e} \n")
