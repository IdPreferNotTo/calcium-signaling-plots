import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import os

import styles as st
import default_parameters as df
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    w = 3.2 * 0.75 * 3
    h = 2.4 * 0.75 * 2
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(w, h), constrained_layout=True)
    for ax in axs:
        ax.remove()
    gridspec = axs[0].get_subplotspec().get_gridspec()
    subfigs = [fig.add_subfigure(gs) for gs in gridspec]

    for row, subfig in enumerate(subfigs):
        if row==0:
            subfig.suptitle(f'mean-driven')
        else:
            subfig.suptitle(f'excitable')

        # create 1x3 subplots per subfig
        axs = subfig.subplots(nrows=1, ncols=3)
        if row==0:
            ax1 = axs[0]
            ax2 = axs[1]
            ax3  =axs[2]
        else:
            ax4 = axs[0]
            ax5 = axs[1]
            ax6  =axs[2]

    ax1.text(0.8, 0.8, r"A$_1$", fontsize=11, transform=ax1.transAxes)
    ax2.text(0.8, 0.8, r"A$_2$", fontsize=11, transform=ax2.transAxes)
    ax3.text(0.8, 0.8, r"A$_3$", fontsize=11, transform=ax3.transAxes)
    ax4.text(0.8, 0.8, r"B$_1$", fontsize=11, transform=ax4.transAxes)
    ax5.text(0.8, 0.8, r"B$_2$", fontsize=11, transform=ax5.transAxes)
    ax6.text(0.8, 0.8, r"B$_3$", fontsize=11, transform=ax6.transAxes)

    taumax = 1000
    taumin = 10
    epsmax = 1.00
    epsmin = 0.01
    size = 41
    N = 10
    n = 5
    m = 3
    home = os.path.expanduser("~")
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]

    cmap_YlGnBu = plt.get_cmap("YlGnBu", 11)
    cmap_YlGnBu_r = plt.get_cmap("YlGnBu_r", 11)
    data = np.loadtxt(home + "/Desktop/out_tau5.0_p0.015.dat")
    tau_ers, eps_ers, ntrs, dTs, rhos =  np.transpose(data)
    tau_ers = tau_ers.reshape(size, size)
    eps_ers = eps_ers.reshape(size, size)
    ntrs = ntrs.reshape(size, size)
    dTs = dTs.reshape(size, size)
    rhos = rhos.reshape(size, size)

    ax1.set_ylabel(r"$\Delta T$ / s")
    ax1.set_xlabel(r"$n_{\rm tr}$")
    ax1.set_ylim([0, 1000])
    ax1.set_xlim([0, 10])
    cm_eps = ax1.pcolormesh(ntrs, dTs, eps_ers,  shading="gouraud", linewidth=0, norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=0.01, vmax=1), rasterized=True, cmap=cmap_YlGnBu)
    cb_cm_eps = fig.colorbar(cm_eps, ax=ax1)
    #ax1.set_title(r"$\varepsilon$", loc="center")
    cb_cm_eps.ax.set_title(r"$\varepsilon$", loc="center")

    ax2.set_ylabel(r"$\Delta T$ / s")
    ax2.set_xlabel(r"$n_{\rm tr}$")
    ax2.set_ylim([0, 1000])
    ax2.set_xlim([0, 10])
    cm_tauers = ax2.pcolormesh(ntrs, dTs, tau_ers,  shading="gouraud", linewidth=0, vmin=0, vmax=1000, rasterized=True, cmap=cmap_YlGnBu)
    cb_cm_tauers = fig.colorbar(cm_tauers, ax=ax2)
    cb_cm_tauers.set_ticks([0, 500, 1000])
    #ax2.set_title(r"$\tau_{\rm er}$ / s", loc="center")
    cb_cm_tauers.ax.set_title(r"$\tau_{\rm er}$ / s", loc="center")

    ax3.set_ylabel(r"$\Delta T$ / s")
    ax3.set_xlabel(r"$n_{\rm tr}$")
    ax3.set_ylim([0, 1000])
    ax3.set_xlim([0, 10])
    cm_rhos = ax3.pcolormesh(ntrs, dTs, rhos,  shading="gouraud", linewidth=0, vmin=-0.5, vmax=0, rasterized=True, cmap=cmap_YlGnBu_r)
    cb_cm_rhos = fig.colorbar(cm_rhos, ax=ax3)
    #ax3.set_title(r"$\rho_1$", loc="center")
    cb_cm_rhos.ax.set_title(r"$\rho_1$", loc="center")

    #cmap_YlGnBu = plt.get_cmap("YlGn", 11)
    #cmap_YlGnBu_r = plt.get_cmap("YlGn_r", 11)
    data = np.loadtxt(home + "/Desktop/out_tau1.0_p0.060.dat")
    tau_ers, eps_ers, ntrs, dTs, rhos =  np.transpose(data)
    size1 = 41
    size2 = 61
    tau_ers = tau_ers.reshape(size1, size2)
    eps_ers = eps_ers.reshape(size1, size2)
    ntrs = ntrs.reshape(size1, size2)
    dTs = dTs.reshape(size1, size2)
    rhos = rhos.reshape(size1, size2)
    ax4.set_ylabel(r"$\Delta T$ / s")
    ax4.set_xlabel(r"$n_{\rm tr}$")
    ax4.set_ylim([0, 1000])
    ax4.set_xlim([0, 5])
    cm_eps = ax4.pcolormesh(ntrs, dTs, eps_ers,  shading="gouraud", linewidth=0, norm=mcolors.LogNorm(vmin=0.001, vmax=1), rasterized=True, cmap=cmap_YlGnBu)
    cb_cm_eps = fig.colorbar(cm_eps, ax=ax4)
    ax4.set_title(r"$\varepsilon$", loc="center")

    ax5.set_ylabel(r"$\Delta T$ / s")
    ax5.set_xlabel(r"$n_{\rm tr}$")
    ax5.set_ylim([0, 1000])
    ax5.set_xlim([0, 5])
    cm_tauers = ax5.pcolormesh(ntrs[0:41, 0:61], dTs[0:41, 0:61], tau_ers[0:41, 0:61], shading="gouraud", linewidth=0, vmin=0, vmax=1000, rasterized=True, cmap=cmap_YlGnBu)
    cb_cm_tauers = fig.colorbar(cm_tauers, ax=ax5)
    cb_cm_tauers.set_ticks([0, 500, 1000])
    ax5.set_title(r"$\tau_{\rm er}$ / s", loc="center")

    ax6.set_ylabel(r"$\Delta T$ / s")
    ax6.set_xlabel(r"$n_{\rm tr}$")
    ax6.set_ylim([0, 1000])
    ax6.set_xlim([0, 5])
    cm_rhos = ax6.pcolormesh(ntrs, dTs, rhos,  shading="gouraud", linewidth=0, vmin=-0.5, vmax=0, rasterized=True, cmap=cmap_YlGnBu_r)
    cb_cm_rhos = fig.colorbar(cm_rhos, ax=ax6)
    ax6.set_title(r"$\rho_1$", loc="center")

    statistics = np.asarray(
        [[40, 0.12, 1.1, 343],
         [24, 0.11, 4.6, 121],
         [21, 0.12, 7.4, 224],
         [13, 0.18, 2.1, 219],
         [25, 0.20, 3.5, 187],
         [19, 0.12, 2.4, 284],
         [29, 0.31, 0.9, 187],
         [29, 0.17, 2.4, 311],
         [29, 0.17, 1.6, 297],
         [29, 0.17, 6.4, 146],
         [21, 0.15, 3.8, 136],
         [22, 0.12, 4.9, 164],
         [16, 0.27, 1.8, 123],
         [35, 0.11, 1.6, 91],
         [20, 0.26, 1.8, 165],
         [21, 0.20, 3.8, 403],
         [41, 0.1, 2.6, 143],
         [46, 0.09, 4.3, 194],
         [46, 0.09, 4.2, 194],
         [31, 0.14, 1.9, 174],
         [36, 0.10, 7.4, 232],
         [36, 0.14, 4.7, 109],
         [15, 0.15, 5.8, 175]])

    t0s = statistics[:, 0]
    cvs = statistics[:, 1]
    ntrs = statistics[:, 2]
    t8s = statistics[:, 3]

    #plt.savefig(home +"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig10.pdf")
    plt.show()