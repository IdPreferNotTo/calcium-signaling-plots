import numpy as np
from scipy.optimize import curve_fit
import os

import functions as fc

if __name__ == "__main__":
    home = os.path.expanduser("~")
    with open(home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/HEK_transient_stationary_statistics.dat", "w") as outfile:
        outfile.write("# idx | T0 | T8 | nTr | std T8 | var T0 | var T8 | var nTr \n")
        with open(home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/spike_times.dat", "r") as infile:
            for idx, line in enumerate(infile):
                if idx not in [7, 10, 17, 19, 21, 22, 24, 25, 26]:
                    data_isi = line[:-2].split(" ")
                    idxs = np.arange(len(data_isi))
                    popt, pcov = curve_fit(fc.exponential_Ti, idxs, data_isi, p0=(100, 150, 2))
                    data_isi_stationary = data_isi[int(1.5 * popt[2] +1):]
                    data_isi_stationary = [float(isi) for isi in data_isi_stationary]
                    variance = np.std(data_isi_stationary)
                    outfile.write(f"{idx+1:d} {popt[0]:.2f} {popt[1]:.2f} {popt[2]:.2f} {variance:.2f} {pcov[0, 0]:.2f} {pcov[1, 1]:.2f} {pcov[2, 2]:.2f}\n")

