import numpy as np
from scipy.optimize import curve_fit
import os

def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))

if __name__ == "__main__":

    tau = 10.5
    j = 0.0146
    N = 10
    home = os.path.expanduser("~")
    folder = home + "/Data/calcium_spikes_markov/Data_adap"
    with open(home + f"/Data/Calcium/data/markov_adap_parameters_tau{tau:.2e}_j{j:.2e}.dat", "w") as outfile:
        outfile.write("# taua | ampa | T0 | T_infty | tau_eff | n \n")
        for i in range(50):
            for k in range(50):
                taua = np.logspace(1, 3, 50)[i]  # 300
                ampa = np.logspace(-2, 1, 50)[k]  # 0.1
                file =  f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N{N:d}_0.dat"
                ISIs = np.loadtxt(folder + file)

                index_ISIs = np.arange(len(ISIs))
                popt, pcov = curve_fit(transient_func, index_ISIs, ISIs, p0=(100, 150, 2))
                outfile.write(f"{taua:.2e} {ampa:.2e} {popt[0]:.2e} {popt[1]:.2e} {popt[2]:.2e} \n")