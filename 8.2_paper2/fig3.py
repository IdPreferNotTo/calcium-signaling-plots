import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc

def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3*2.25, 4))
    gs = gridspec.GridSpec(1, 2)
    gs1 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[0], hspace=1.5)
    gs2 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[1], hspace=1.5)
    ax0 = fig.add_subplot(gs1[0:2])
    ax1 = fig.add_subplot(gs1[2])
    ax2 = fig.add_subplot(gs1[3:5])

    ax3 = fig.add_subplot(gs2[0:2])
    ax4 = fig.add_subplot(gs2[2])
    ax5 = fig.add_subplot(gs2[3:5])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3, ax4, ax5])

    ax0.text(0.1, 0.95, "A$_{i}$", fontsize=13, transform=ax0.transAxes, va='top')
    ax1.text(0.1, 0.95, "A$_{ii}$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, "A$_{iii}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.1, 0.95, "B$_{i}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.1, 0.95, "B$_{ii}$", fontsize=13, transform=ax4.transAxes, va='top')
    ax5.text(0.1, 0.95, "B$_{iii}$", fontsize=13, transform=ax5.transAxes, va='top')

    taua = 244
    ampa = 0.11
    home = os.path.expanduser("~")
    isis = []
    for i in range(10):
        folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/P2F3"
        file_isi = f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_{i:d}.dat"
        isis_full = np.loadtxt(folder + file_isi)
        for I in isis_full[100:]:
            isis.append(I)
    mean = np.mean(isis)
    std = np.std(isis)
    var= np.var(isis)
    cv = std/mean
    print(mean)

    ks = np.arange(1, 5)
    rho_ks = [[] for _ in ks]
    for I in fc.chunks(isis, 100):
        for k in ks:
            var_k_I = fc.k_corr(I, I, k)
            var_I = np.var(I)
            rho_k = var_k_I/var_I
            rho_ks[k-1].append(rho_k)

    mean_rho_ks = [np.mean(rho_k) for rho_k in rho_ks]
    sum_rho_ks = sum(mean_rho_ks)

    ws = np.logspace(-3, 0, 100)
    spectrum = fc.power_spectrum_isis(ws, isis, 2000)
    isis_shuffled = list(isis)
    np.random.shuffle(isis_shuffled)
    spectrum_shuffle = fc.power_spectrum_isis(ws, isis_shuffled, 2000)

    nr_isis = 40
    index_isis = np.arange(nr_isis)
    popt, pcov = curve_fit(transient_func, index_isis, isis_full[0:nr_isis], p0=(100, 150, 2))
    ax0.set_xlim([0, nr_isis])
    ax0.set_xlabel("$i$")
    ax0.set_ylim([0, 1.5 * popt[1]])
    ax0.set_ylabel("$T_i$ / s")
    ax0.scatter(np.arange(0, nr_isis), isis_full[0:nr_isis], fc="w", ec=st.colors[1], s=20, zorder=3)
    isis_fit = popt[0] * np.exp(-index_isis / popt[2]) + popt[1] * (1. - np.exp(-index_isis / popt[2]))
    ax0.plot(index_isis, isis_fit, lw=1, c="k")
    ax0.axhline(popt[0], ls=":", lw=1, c="k")
    ax0.axhline(popt[1], ls=":", lw=1, c="k")
    ax0.text(nr_isis/2, popt[0]*1.2, "$T_0$", ha="center")
    ax0.text(nr_isis/2, popt[1]*1.2, "$T_\infty$", ha="center")
    ax0.fill_between([0, popt[2]], [0, 0], [200, 200], color="C7", alpha=0.5, zorder=1)

    ax1.set_ylim([-0.3, 0.05])
    ax1.set_xlabel("$k$")
    ax1.set_ylabel(r"$\rho_k$")
    ax1.scatter(ks, [np.mean(rho_k) for rho_k in rho_ks], fc="w", ec=st.colors[1], s=20, zorder=3)
    for k, rho_k in enumerate(rho_ks):
        mean_rho_k = np.mean(rho_k)
        std_error = np.std(rho_k)/np.sqrt(len(rho_k))
        ax1.plot([k+1, k+1], [mean_rho_k - std_error, mean_rho_k + std_error], c=st.colors[1])
    ax1.axhline(0, ls=":", c="C7")

    ax2.set_xlabel("$\omega$")
    ax2.set_ylabel("$S(\omega)$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim([0.001, 2])
    ax2.plot(ws, spectrum, color=st.colors[1], label="original")
    ax2.plot(ws, spectrum_shuffle, color=st.colors[0], label="shuffled")
    ax2.legend(frameon=False, loc=4)

    ax2.axhline((1. / mean), ls="--", c="C7")
    ax2.axhline((1. / mean) * cv ** 2, ls="--", c="C7")
    ax2.axhline((1. / mean) * cv ** 2 * (1 + 2*sum_rho_ks), ls="--", c="C7")

    #------------------------------------------------------------------------------------------------------------------

    taua = 829
    ampa = 0.0309
    home = os.path.expanduser("~")
    isis = []
    for i in range(10):
        folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/P2F3"
        file_isi = f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_{i:d}.dat"
        isis_full = np.loadtxt(folder + file_isi)
        for I in isis_full[100:]:
            isis.append(I)
    mean = np.mean(isis)
    std = np.std(isis)
    var= np.var(isis)
    cv = std/mean
    print(mean)

    ks = np.arange(1, 5)
    rho_ks = [[] for _ in ks]
    for I in fc.chunks(isis, 100):
        for k in ks:
            var_k_I = fc.k_corr(I, I, k)
            var_I = np.var(I)
            rho_k = var_k_I/var_I
            rho_ks[k-1].append(rho_k)

    mean_rho_ks = [np.mean(rho_k) for rho_k in rho_ks]
    sum_rho_ks = sum(mean_rho_ks)

    ws = np.logspace(-3, 0, 100)
    spectrum = fc.power_spectrum_isis(ws, isis, 2000)
    isis_shuffled = list(isis)
    np.random.shuffle(isis_shuffled)
    spectrum_shuffle = fc.power_spectrum_isis(ws, isis_shuffled, 2000)

    nr_isis = 40
    index_isis = np.arange(nr_isis)
    popt, pcov = curve_fit(transient_func, index_isis, isis_full[0:nr_isis], p0=(100, 150, 2))
    ax3.set_xlim([0, nr_isis])
    ax3.set_xlabel("$i$")
    ax3.set_ylim([0, 1.5 * popt[1]])
    ax3.set_ylabel("$T_i$ / s")
    ax3.scatter(np.arange(0, nr_isis), isis_full[0:nr_isis], fc="w", ec=st.colors[1], s=20, zorder=3)
    isis_fit = popt[0] * np.exp(-index_isis / popt[2]) + popt[1] * (1. - np.exp(-index_isis / popt[2]))
    ax3.plot(index_isis, isis_fit, lw=1, c="k")
    ax3.axhline(popt[0], ls=":", lw=1, c="k")
    ax3.axhline(popt[1], ls=":", lw=1, c="k")
    ax3.text(nr_isis/2, popt[0]*1.2, "$T_0$", ha="center")
    ax3.text(nr_isis/2, popt[1]*1.2, "$T_\infty$", ha="center")
    ax3.fill_between([0, popt[2]], [0, 0], [200, 200], color="C7", alpha=0.5, zorder=1)

    ax4.set_ylim([-0.3, 0.05])
    ax4.set_xlabel("$k$")
    ax4.set_ylabel(r"$\rho_k$")
    ax4.scatter(ks, [np.mean(rho_k) for rho_k in rho_ks], fc="w", ec=st.colors[1], s=20, zorder=3)
    for k, rho_k in enumerate(rho_ks):
        mean_rho_k = np.mean(rho_k)
        std_error = np.std(rho_k)/np.sqrt(len(rho_k))
        ax4.plot([k+1, k+1], [mean_rho_k - std_error, mean_rho_k + std_error], c=st.colors[1])
    ax4.axhline(0, ls=":", c="C7")

    ax5.set_xlabel("$\omega$")
    ax5.set_ylabel("$S(\omega)$")
    ax5.set_xscale("log")
    ax5.set_yscale("log")
    ax5.set_xlim([0.001, 2])
    ax5.plot(ws, spectrum, color=st.colors[1], label="original")
    ax5.plot(ws, spectrum_shuffle, color=st.colors[0], label="shuffled")
    ax5.legend(frameon=False, loc=4)

    ax5.axhline((1. / mean), ls="--", c="C7")
    ax5.axhline((1. / mean) * cv ** 2, ls="--", c="C7")
    ax5.axhline((1. / mean) * cv ** 2 * (1 + 2*sum_rho_ks), ls="--", c="C7")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig3.png", transparent=True)
    plt.show()