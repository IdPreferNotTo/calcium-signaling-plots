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
    file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
    data = np.loadtxt(file_str)


    rho_1s = []
    i_effs = []
    n = len(data[0])
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

        ISIs = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]
        idx_ISIs = [i for i in range(len(ISIs))]
        popt, pcov = curve_fit(transient_func, idx_ISIs, ISIs, p0=(ISIs[0], ISIs[-1], 10), bounds=[[0, 0, 0], [np.inf, np.inf, np.inf]])

        i_stat = int(2*popt[2])
        if popt[2] < 15 and j not in [7, 10, 17, 19, 21, 22, 24, 25, 26] and len(ISIs[i_stat:]) > 10:

            var_0 = fc.k_corr(ISIs[i_stat:], ISIs[i_stat:], 0)
            var_1 = fc.k_corr(ISIs[i_stat:], ISIs[i_stat:], 1)
            rho_1 = var_1/var_0
            rho_1s.append(rho_1)
            i_effs.append(popt[2])
            print(j, rho_1)

        st.set_default_plot_style()
        fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])
        st.remove_top_right_axis([ax1])

        Tis = [popt[0] * np.exp(-idx/popt[2]) + popt[1] * (1 - np.exp(-idx/popt[2])) for idx in idx_ISIs]
        ax1.plot(idx_ISIs, Tis, c="k", lw=1, label=f"$T_0$ = {popt[0]:.0f}" + "\n" + f"$T_\infty$ = {popt[1]:.0f}" + "\n" + rf"$\iota_{{\rm eff}}$ = {popt[2]:.1f}")
        ax1.set_ylim([0, 1.1*max(ISIs)])
        ax1.set_xlim([0, len(ISIs)])
        ax1.set_xlabel("$i$")
        ax1.set_ylabel(r"$T_i$ / s")
        ax1.scatter(idx_ISIs, ISIs, s=20, fc="w", ec=st.colors[4])
        ax1.axvspan(0, 2*popt[2], facecolor="C7", alpha=0.5, zorder=0)
        ax1.axhline(popt[0], ls=":", c="C7")
        ax1.axhline(popt[1], ls=":", c="C7")
        legend = ax1.legend(loc=4, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 10})
        legend.get_frame().set_linewidth(0.5)
        plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/Adaptation/HEK2_bapta_{:d}_adaptation.png".format(j))
        #plt.show()
        plt.close()

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax1])
    ax1.scatter(i_effs, rho_1s, s=20, fc="w", ec=st.colors[4])
    ax1.set_xlabel(r"$\iota_{\rm eff}$")
    ax1.set_ylabel(r"$\rho_1$")
    plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/Adaptation/HEK2_bapta_iota_rho.png".format(j))
    plt.show()


