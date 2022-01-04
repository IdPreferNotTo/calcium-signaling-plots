import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os
from scipy.optimize import curve_fit

import functions as fc

def lin_func(x, m):
    return m*x

if __name__ == "__main__":
    home = os.path.expanduser("~")
    file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
    data = np.loadtxt(file_str)
    n = len(data[0])
    for j in range(1, n):
        ts = [x[0] for x in data]
        cas = [x[j] for x in data]
        spiking: bool = False
        t_tmp: float = 0
        spike_times = []
        if j == 2:
            for t, ca in zip(ts, cas):
                if t < 3700:
                    if ca > 0.35 and not spiking:
                        spike_times.append(t)
                        spiking = True
                    if ca < 0.35 and spiking:
                        spiking = False
        else:
            for t, ca in zip(ts, cas):
                if t < 3700:
                    if ca > 0.4 and not spiking:
                        spike_times.append(t)
                        spiking = True
                    if ca < 0.4 and spiking:
                        spiking = False
        print(j)
        ISIs = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]
        t_ISIs = []
        T = 0
        for isi in ISIs:
            T += isi
            t_ISIs.append(T)

        fc.set_default_plot_style()
        fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
        gs = gridspec.GridSpec(2, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[0, 1])
        ax4 = fig.add_subplot(gs[1, 1])
        fc.remove_top_right_axis([ax1, ax2, ax3, ax4])

        ax1.set_xlabel("$i$")
        ax1.set_ylabel("$T_\infty - T_i$")

        ax3.set_xlabel("$t$")
        ax3.set_ylabel("$T_\infty - T_i(t)$")

        ax2.set_xlabel("$T_\infty - T_{i}$")
        ax2.set_ylabel("$T_\infty - T_{i+1}$")

        ax4.set_xlabel("$i$")
        ax4.set_ylabel("$T_{i+1} - T_i$")

        ISIs_ax1 = [ISIs[-1] - I for I in ISIs]
        ax1.scatter(range(len(ISIs_ax1)), ISIs_ax1, s=10, fc="w", ec="C0")

        ax3.scatter(t_ISIs, ISIs_ax1, s=10, fc="w", ec="C0")

        ax2.scatter(ISIs_ax1[0:-1], ISIs_ax1[1:], s=10, fc="w", ec="C0")

        ISIs_ax4 =  [I2 - I1 for I1, I2 in zip(ISIs[0:-1], ISIs[1:])]
        ax4.scatter(range(len(ISIs_ax4)), ISIs_ax4, s=10, fc="w", ec="C0")
        plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/Adaptation/HEK2_bapta_{:d}_stochastic_map.png".format(j))
        plt.show()
