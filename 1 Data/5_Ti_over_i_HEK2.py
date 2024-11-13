import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from scipy.optimize import curve_fit

import functions as fc
import styles as st

def transient_func(i, T0, T8, iota):
    return T0*np.exp(-i/iota) + T8*(1. - np.exp(-i/iota))

if __name__ == "__main__":
    thresholds = [0.4, 0.35, 0.4, 0.35, 0.4, 0.35, 0.35, 0.35, 0.35, 0.4, 0.45, 0.35, 0.35, 0.45, 0.4, 0.35, 0.4, 0.45,
                  0.35, 0.35, 0.5, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.4, 0.45, 0.45, 0.45,
                  0.35, 0.4, 0.4]

    home = os.path.expanduser("~")
    file_str = home + "/Data/calcium/experimental/Spikes/HEK/HEK2_bapta_ratio.dat"
    data = np.loadtxt(file_str)


    rho_1s = []
    i_effs = []
    n = len(data[0])
    with open(home + f"/Data/calcium/experimental/Spikes/HEK/HEK2/spike_times.dat", "w") as f:
        for j in range(1, n):
            ts = [x[0] for x in data]
            cas = [x[j] for x in data]

            spiking: bool = False
            t_tmp: float = 0
            spike_times = []
            for t, ca in zip(ts, cas):
                if t < 3700:
                    if ca > thresholds[j-1] and not spiking:
                        spike_times.append(t)
                        spiking = True
                    if ca < thresholds[j-1] and spiking:
                        spiking = False

            ISIs = []
            t0 = 0
            for t in spike_times:
                ISIs.append(t - t0)
                t0 = t
            with open(home + f"/Data/calcium/experimental/Spikes/HEK/HEK2/spike_times_{j:d}.dat", "w") as f2:
                for ISI in ISIs[:-1]:
                    f.write(f"{ISI:.1f} ")
                    f2.write(f"{ISI:.1f} ")
                f.write(f"{ISIs[-1]:.1f}") # no closing space
                f2.write(f"{ISIs[-1]:.1f}")
                f.write("\n")


            #ISIs = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]

            idx_ISIs = [i for i in range(len(ISIs))]
            popt, pcov = curve_fit(fc.exponential_Ti, idx_ISIs, ISIs, bounds=([0, 0 ,0], [np.inf, np.inf, np.inf]))
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            perr = np.sqrt(np.diag(pcov))

            i_stat = int(2*popt[2])
            cv = np.nan
            if popt[2] < 15 and j not in [10, 17, 21, 22, 24, 26]:
                std_isi = np.std(ISIs[i_stat:])
                mean_isi = np.mean(ISIs[i_stat:])
                cv = std_isi/mean_isi
                print(j, std_isi/mean_isi, std_isi/popt[1])


            st.set_default_plot_style()
            fig = plt.figure(tight_layout=True, figsize=(4, 2.5))
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0])
            st.remove_top_right_axis([ax1])

            Tis = [popt[0] * np.exp(-idx/popt[2]) + popt[1] * (1 - np.exp(-idx/popt[2])) for idx in idx_ISIs]
            ax1.plot(idx_ISIs, Tis, c="k", lw=1, label=f"$T_0$ = {popt[0]:.0f} $\pm$ {perr[0]:.0f}" + "\n" + f"$T_\infty$ = {popt[1]:.0f} $\pm$ {perr[1]:.0f}" + "\n" + rf"$n_{{\rm tr}}$ = {popt[2]:.1f} $\pm$ {perr[2]:.1f}" + "\n" + rf"$CV = {cv:.2f}$")
            ax1.set_ylim([0, 1.1*max(ISIs)])
            ax1.set_xlim([0, len(ISIs)])
            ax1.set_xlabel("$i$")
            ax1.set_ylabel(r"$T_i$ / s")
            ax1.scatter(idx_ISIs, ISIs, s=20, fc="w", ec=st.colors[7])
            ax1.axvspan(0, popt[2], facecolor="C7", alpha=0.5, zorder=0)
            ax1.axhline(popt[0], ls=":", c="C7")
            ax1.axhline(popt[1], ls=":", c="C7")
            legend = ax1.legend(loc=4, fancybox=False, edgecolor="k", framealpha=0.0, prop={'size': 10})
            legend.get_frame().set_linewidth(0.5)

            plt.savefig(home + "/Data/calcium/experimental/Spikes/HEK/Plots2/Adaptation/HEK2_bapta_{:d}_adaptationpng".format(j))
            plt.show()
            plt.close()
