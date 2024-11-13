import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.66 * 1.25
    fig, axs = plt.subplots(2, 2, layout="constrained", figsize=(w, h))
    ax1 = axs[1,1]
    ax2 = axs[1,0]
    ax3 = axs[0,1]
    ax4 = axs[0,0]

    ax1.text(0.05, 0.95, r"B$_2$", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B$_1$", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"A$_2$", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"A$_1$", fontsize=11, transform=ax4.transAxes, va='top')

    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    ax1.set_ylabel(r"$\rho_1$")
    ax1.set_xlabel(r"$n_{\rm tr}$")
    ax2.set_ylabel(r"$\rho_1$")
    ax2.set_xlabel(r"$\Delta T$ / s")
    ax3.set_ylabel(r"$\rho_1$")
    ax3.set_xlabel(r"$n_{\rm tr}$")
    ax4.set_ylabel(r"$\rho_1$")
    ax4.set_xlabel(r"$\Delta T$ / s")
    ax1.set_xlim([0, 5.5])
    ax3.set_xlim([0, 5.5])
    ax2.set_xlim([0, 800])
    ax4.set_xlim([0, 800])

    for ax_row in axs:
        for ax in ax_row:
            ax.set_ylim([-0.5, 0.15])
            ax.axhline(0., ls=":", c="k")

    # Parameters
    tau = 5.0
    p = 0.015
    eps_er_fix = 0.03
    tau_er_fix = 300
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)
    cmap_YlGnBu = plt.get_cmap("YlGnBu", 21)

    n_trs = []
    dTs = []
    p1s = []

    n_tr_crit = 0
    dT_crit = 0
    for eps_er in eps_ers:
        data_isi_stationary = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er_fix, ampa=eps_er)
        var = fc.k_corr(data_isi_stationary, data_isi_stationary, k=0)
        p1 = fc.k_corr(data_isi_stationary, data_isi_stationary, k=1)
        p2 = fc.k_corr(data_isi_stationary, data_isi_stationary, k=2)
        p1s.append(p1/var)

        data_isi_transient = df.load_spike_times_markov_transient(tau, p, tau_er_fix, eps_er)
        rows, cols = data_isi_transient.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        means_Tidx = [np.mean(data_isi_transient[:, idx]) for idx in idxs]  # calculate the mean column wise
        popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        T0 = popt[0]
        T8 = popt[1]
        n_tr = popt[2]
        n_trs.append(n_tr)
        dTs.append(T8 - T0)
        if eps_er == eps_ers[8]:
            print(eps_er)
            n_tr_crit =  n_tr
            dT_crit = T8 - T0

    ax1.axvspan(n_tr_crit, 6, alpha=0.3, color="C7")
    ax2.axvspan(0, dT_crit, alpha=0.3, color="C7")
    ax1.set_xlim([0, 5.5])
    ax2.set_xlim([0, 800])
    #ax1.axvline(n_tr_crit, ls="--", c="k")
    #ax2.axvline(dT_crit, ls="--", c="k")
    c = plt.cm.YlGnBu(eps_ers)
    scat1 = ax1.scatter(n_trs, p1s, c=eps_ers, s=20, norm=mcolors.LogNorm(vmin=0.01, vmax=1.), cmap=cmap_YlGnBu, zorder=3)
    scat2 = ax2.scatter(dTs, p1s, c=eps_ers, s=20, norm=mcolors.LogNorm(vmin=0.01, vmax=1.), cmap=cmap_YlGnBu, zorder=3)
    scat1.set_edgecolor(st.colors[0])
    scat2.set_edgecolor(st.colors[0])
    cbar = plt.colorbar(scat1, ax=[ax1, ax2])
    cbar.set_label(r"$\varepsilon$")

    n_trs = []
    dTs = []
    p1s = []
    for tau_er in tau_ers:
        data_isi_stationary = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er_fix)
        var = fc.k_corr(data_isi_stationary, data_isi_stationary, k=0)
        p1 = fc.k_corr(data_isi_stationary, data_isi_stationary, k=1)
        p2 = fc.k_corr(data_isi_stationary, data_isi_stationary, k=2)
        p1s.append(p1/var)

        data_isi_transient = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er_fix)
        rows, cols = data_isi_transient.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        means_Tidx = [np.mean(data_isi_transient[:, idx]) for idx in idxs]  # calculate the mean column wise
        popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=(100, 150, 2))
        T0 = popt[0]
        T8 = popt[1]
        n_tr = popt[2]
        n_trs.append(n_tr)
        dTs.append(T8 - T0)

    c = plt.cm.YlGnBu(tau_ers)
    scat3 = ax3.scatter(n_trs, p1s, c=tau_ers, s=20, norm=mcolors.LogNorm(vmin=10., vmax=1000), cmap=cmap_YlGnBu, zorder=3)
    scat4 = ax4.scatter(dTs, p1s, c=tau_ers, s=20, norm=mcolors.LogNorm(vmin=10., vmax=1000), cmap=cmap_YlGnBu, zorder=3)
    scat3.set_edgecolor(st.colors[0])
    scat4.set_edgecolor(st.colors[0])
    cbar = plt.colorbar(scat3, ax=[ax3, ax4])
    cbar.set_label(r"$\tau_{\rm er}$ / s")

home = os.path.expanduser("~")
plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig9.pdf", dpi=300, transparent=True)
plt.show()
