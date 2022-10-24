import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
import os

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    home = os.path.expanduser("~")
    taus = [5, 1]
    ps = [0.015, 0.06]
    IP3s = [1.]
    for tau, p in zip(taus, ps):
        print(tau, p)
        for IP3 in IP3s:
            isis_markov = df.load_spike_times_markov(tau, p, cer=False, ip3=IP3)
            with open(home + f"/Data/calcium_spikes_theory/Fano_factor_ip{IP3:.2f}_tau{tau:.2e}_j{p:.2e}.dat",
                      "w") as file:
                file.write("#Delta T mean-ISI std-ISI Fano-ISI \n")
                dts = np.logspace(-1, 3, 1000)
                mean_counts = []
                var_counts = []
                for dt in dts:
                    print(dt)
                    mean_count, var_count = fc.fano_factor_interspike_intervals(isis_markov, dt)
                    file.write(f"{dt:.2e} {mean_count:.2e} {var_count:.2e} {var_count/mean_count:.2e} \n")
