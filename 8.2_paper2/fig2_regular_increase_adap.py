import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

import styles as st
import functions as fc

def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))

def inverse_gaussian(T, CV):
    ps = []
    ts = np.linspace(0, 2*T, 500)
    for t in ts[1:]:
        p = np.sqrt(T / (2 * np.pi * (CV**2) * (t ** 3))) * np.exp(-(t - T) ** 2 / (2 * T * (CV**2) * t))
        ps.append(p)
    return ts[1:], ps

if __name__ == "__main__":
    ep_s = [0.05, 0.1, 0.2]
    taua_s = [100, 200, 300]


    home = os.path.expanduser("~")

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(9, 6))
    gs = gridspec.GridSpec(nrows=3, ncols=3)

    pattern = "regular"
    if pattern == "regular":
        tau = 10.  # 10.5
        p = 0.015  # 0.0146
    elif pattern == "medium":
        tau = 1.7
        p = 0.06
    else:
        tau = 0.5
        p = 0.2
    axiss = []
    for j, e in enumerate(ep_s):
        axis = []
        for i , taua in enumerate(taua_s):
            ax = fig.add_subplot(gs[i, j])
            axis.append(ax)
            if j == 2:
                axins = inset_axes(ax, width=1.0, height=0.5, loc=4)
            else:
                axins = inset_axes(ax, width=1.0, height=0.5, loc=1)
            st.remove_top_right_axis([ax])
            if pattern == "irregular":
                if i == 0:
                    ax.set_ylim([0, 150]) #ax.set_ylim([0, 100])
                elif i == 1:
                    ax.set_ylim([0, 200]) #ax.set_ylim([0, 150])
                else:
                    ax.set_ylim([0, 300]) #ax.set_ylim([0, 200])
            else:
                if i == 0:
                    ax.set_ylim([0, 100]) #ax.set_ylim([0, 100])
                elif i == 1:
                    ax.set_ylim([0, 150]) #ax.set_ylim([0, 150])
                else:
                    ax.set_ylim([0, 200]) #ax.set_ylim([0, 200])
            folder = home + "/Data/calcium_spikes_markov/Data_adap"
            file_isi = f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{e:.2e}_tau{tau:.2e}_j{p:.2e}_N10_0.dat"
            isis = np.loadtxt(folder + file_isi)
            rho_1 = fc.k_corr(isis[100:], isis[100:], 1)/np.var(isis[100:])

            axins.hist(isis[20:], bins=20, color=st.colors[1], density=True, alpha=0.6)
            axins.set_xlabel("")
            axins.set_xticks([])
            axins.set_ylabel("$p(T)$")
            axins.set_yticks([])
            axins.set_title(rf"$\rho_1 = {rho_1:.2f}$")
            meanT = np.mean(isis[20:])
            cvT = np.std(isis[20:])/meanT
            ts_inv_gaussian, ps_inv_gaussian = inverse_gaussian(meanT, cvT)
            axins.plot(ts_inv_gaussian, ps_inv_gaussian, c="k")

            nr_isis = 20
            ax.set_xlim([0, nr_isis])
            # ax_2.set_yticks([0, 50, 100, 150])
            if j == 0:
                ax.set_ylabel("$T_i$ / s")
            if i ==2:
                ax.set_xlabel("$i$")
            ax.scatter(np.arange(0, nr_isis), isis[0:nr_isis], fc="w", ec=st.colors[1], s=20, zorder=3)
            index_isis = np.arange(len(isis))
            popt, pcov = curve_fit(transient_func, index_isis, isis, p0=(100, 150, 2))
            ax.axvspan(0, popt[2], alpha=0.3, color="C7")

            #print(popt)
            ISI_fit = popt[0] * np.exp(-index_isis / popt[2]) + popt[1] * (1. - np.exp(-index_isis / popt[2]))
            ax.plot(index_isis, ISI_fit, lw=1, c="k")
            ax.axhline(popt[0], ls=":", lw=1, c="k")
            ax.axhline(popt[1], ls=":", lw=1, c="k")
            #ax.text(nr_isis / 2, popt[0] * 1.5, "$T_0$", ha="center")
            #ax.text(nr_isis / 2, popt[1] * 1.2, "$T_\infty$", ha="center")


        axiss.append(axis)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig2_{pattern:s}.pdf", transparent=True)
    plt.show()