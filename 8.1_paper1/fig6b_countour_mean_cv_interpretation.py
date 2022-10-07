import numpy as np
from matplotlib import gridspec, cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import styles as st

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(figsize=(8, 4))
    gs = gridspec.GridSpec(2, 3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])
    axes_mean = [ax1, ax2, ax3]
    axes_cv = [ax4, ax5, ax6]

    ax1.text(0.1, 0.1, r"A$_1$", fontsize=13, transform=ax1.transAxes)
    ax2.text(0.1, 0.1, r"A$_2$", fontsize=13, transform=ax2.transAxes)
    ax3.text(0.1, 0.1, r"A$_3$", fontsize=13, transform=ax3.transAxes)
    ax4.text(0.1, 0.1, r"B$_1$", fontsize=13, transform=ax4.transAxes)
    ax5.text(0.1, 0.1, r"B$_2$", fontsize=13, transform=ax5.transAxes)
    ax6.text(0.1, 0.1, r"B$_3$", fontsize=13, transform=ax6.transAxes)
    colors = st.Colors()
    cmap_coolwarm = plt.get_cmap("coolwarm", 10)

    size = 41
    N = 10
    n = 5
    m = 3
    home = os.path.expanduser("~")
    data_markov = np.loadtxt(home + f"/Data/calcium_spikes_theory/markov_ca_mean_CV_K{N:d}_N{n:d}_no_adap.dat")
    taus_m, js_m, Ts_m, cvs_m, num_m = np.transpose(data_markov)

    taus = np.logspace(0, 2, size)
    js = np.logspace(-3, -1, size)

    r_cls = 50
    r_ref = 20
    r_opn_single = 0.1
    ropnmax = 5 * r_opn_single * ((1. + np.power(0.20, 3)) / np.power(0.20, 3))
    r_opn_ct = ropnmax * (np.power(0.5, 3) / (1. + np.power(0.5, 3)))
    mean_puff = (6) * (7) / (6 * r_cls)
    tau_tot = 1 / r_opn_ct + 2 / r_ref + 6 / (2 * r_cls)
    mean_puff_ct = mean_puff / tau_tot
    ps_bif = []
    for tau in taus:
        ps_bif.append((0.5 - 0.2) / (10 * mean_puff_ct * tau))

    ISI_markov = np.empty([size, size])
    CV_markov = np.empty([size, size])
    for k, T in enumerate(Ts_m):
        if T >= 1_000:
            ISI_markov[k // size, k % size] = np.nan
        else:
            ISI_markov[k // size, k % size] = T
    for k, cv in enumerate(cvs_m):
        if cv == 1:
            CV_markov[k // size, k % size] = np.nan
        else:
            CV_markov[k // size, k % size] = cv

    pmax = 0.1
    pmin = 0.001
    tmax = 100
    tmin = 1.0
    interpretations = ["Ito", "Stratonovich", "Klimontovich-Hänggi"]
    for i, interpretation in enumerate(interpretations):
        ax_mean = axes_mean[i]
        ax_cv = axes_cv[i]
        ax_mean.set_xscale("log")
        ax_mean.set_yscale("log")
        ax_mean.set_xlim([tmin, tmax])
        ax_mean.set_ylim([pmin, pmax])
        ax_mean.set_xlim([tmin, tmax])
        ax_mean.set_ylim([pmin, pmax])

        ax_cv.set_xlabel(r"$\tau$")
        ax_cv.set_xscale("log")
        ax_cv.set_yscale("log")
        ax_cv.set_xlim([tmin, tmax])
        ax_cv.set_ylim([pmin, pmax])
        ax_cv.set_xlim([tmin, tmax])
        ax_cv.set_ylim([pmin, pmax])

        ax_mean.plot(taus, ps_bif, ls="--", c="k")
        ax_cv.plot(taus, ps_bif, ls="--", c="k")
        if interpretation == "Ito":
            interpret = "ito"
            ax_mean.set_ylabel("$p$")
            ax_mean.set_xticklabels([])
            ax_cv.set_ylabel("$p$")
            a = 0.0
        elif interpretation == "Stratonovich":
            interpret = "strat"
            ax_mean.set_xticklabels([])
            ax_mean.set_yticklabels([])
            ax_cv.set_yticklabels([])
            a = 0.5
        elif interpretation == "Klimontovich-Hänggi":
            interpret = "haenggi"
            ax_mean.set_xticklabels([])
            ax_mean.set_yticklabels([])
            ax_cv.set_yticklabels([])
            a = 1.0


        data_langevin = np.loadtxt(home + f"/Data/calcium_spikes_theory/langevin_{interpret:s}_ca_mean_CV_K{N:d}_N{n:d}_no_adap.dat")

        taus_m, js_m, Ts_m, cvs_m, num_m = np.transpose(data_markov)
        taus_l, js_l, Ts_l, cvs_l, num_l = np.transpose(data_langevin)

        ISI_langevin = np.empty([size, size])
        CV_langevin = np.empty([size, size])
        dif_ISI = np.empty([size, size])
        dif_CV = np.empty([size, size])

        for k, T in enumerate(Ts_l):
            if T >= 1_000:
                ISI_langevin[k // size, k % size] = np.nan
                dif_ISI[k // size, k % size] = np.nan
            else:
                ISI_langevin[k // size, k % size] = T
                dif_ISI[k // size, k % size] = -(ISI_markov[k // size, k % size] - T) / ISI_markov[k // size, k % size]
        for k, cv in enumerate(cvs_l):
            if cv == 1:
                CV_langevin[k // size, k % size] = np.nan
                dif_CV[k // size, k % size] = np.nan
            else:
                CV_langevin[k // size, k % size] = cv
                dif_CV[k // size, k % size] = -(CV_markov[k // size, k % size] - cv) / CV_markov[k // size, k % size]

        print(np.nanmin(dif_ISI), np.nanmax(dif_ISI))
        print(np.nanmin(dif_CV), np.nanmax(dif_CV))
        cs1 = ax_mean.pcolormesh(taus, js, dif_ISI, linewidth=0, rasterized=True, shading='gouraud', vmin=-0.5, vmax=0.5,
                                      cmap=cmap_coolwarm)
        cs2 = ax_cv.pcolormesh(taus, js, dif_CV, linewidth=0, rasterized=True, shading='gouraud', vmin=-0.5, vmax=0.5, cmap=cmap_coolwarm)
        ax_mean.set_title(f"{interpretation:s}")

    cbar1 = fig.colorbar(cs1, ax=axes_mean, orientation='vertical')
    cbar1.set_ticks([-0.5, 0., 0.5])
    cbar1.set_label(r"$(\langle \widetilde{T} \rangle - \langle T \rangle )/ \langle T \rangle$", loc="center")

    cbar2 = fig.colorbar(cs2, ax=axes_cv, orientation='vertical')
    cbar2.set_ticks([-0.5, 0., 0.5])
    cbar2.set_label(r"$(\widetilde{CV}_T - CV_T )/ CV_T$", loc="center")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig6b.pdf", transparent=True)
    plt.show()


