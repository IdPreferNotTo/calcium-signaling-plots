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
    fig = plt.figure(constrained_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    ax1.set_ylabel(r"$\rho_1$")
    ax1.set_xlabel(r"$n_{\rm tr}$")
    ax2.set_ylabel(r"$\rho_1$")
    ax2.set_xlabel(r"$\Delta T$")
    ax3.set_ylabel(r"$\rho_1$")
    ax3.set_xlabel(r"$n_{\rm tr}$")
    ax4.set_ylabel(r"$\rho_1$")
    ax4.set_xlabel(r"$\Delta T$")

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
        popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=(100, 150, 2))
        T0 = popt[0]
        T8 = popt[1]
        n_tr = popt[2]
        n_trs.append(n_tr)
        dTs.append(T8 - T0)


    c = plt.cm.YlGnBu(eps_ers)
    scat1 = ax1.scatter(n_trs, p1s, c=eps_ers, s=20, norm=mcolors.LogNorm(vmin=0.01, vmax=1.), cmap=cmap_YlGnBu, zorder=3)
    scat2 = ax2.scatter(dTs, p1s, c=eps_ers, s=20, norm=mcolors.LogNorm(vmin=0.01, vmax=1.), cmap=cmap_YlGnBu, zorder=3)
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
    cbar = plt.colorbar(scat3, ax=[ax3, ax4])
    cbar.set_label(r"$\tau_{\rm er}$")

home = os.path.expanduser("~")
plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig8.pdf", transparent=True)
plt.show()
