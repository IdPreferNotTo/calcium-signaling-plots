import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, cm
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import styles as st

if __name__ == "__main__":

    st.set_default_plot_style()
    fig = plt.figure(figsize=(9., 4.))
    w0 = 0.075
    h0 = 0.15

    w = 0.17
    h = 0.325
    w_cb = w/10
    dw1 = 0.11
    dw2 = 0.02
    dh = 0.1

    ax_mean_markov = fig.add_axes([w0, h0 + h + dh, w, h])
    ax_cv_markov = fig.add_axes([w0, h0, w, h])

    size = 61
    N = 10
    n = 5
    home = os.path.expanduser("~")
    data_markov = np.loadtxt(home + f"/Data/calcium_spikes_theory/markov_ca_mean_CV_K{N:d}_N{n:d}_no_adap.dat")
    taus_m, js_m, Ts_m, cvs_m, num_m = np.transpose(data_markov)

    taus = np.logspace(-1, 2, size)
    js = np.logspace(-3, 0, size)

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
    # st.remove_top_right_axis([ax11, ax12, ax13, ax14, ax21, ax22, ax23, ax24])

    cmap_cividis = plt.get_cmap("cividis", 8)
    cmap_coolwarm = plt.get_cmap("coolwarm", 8)
    cs_mean_markov = ax_mean_markov.pcolormesh(taus, js, ISI_markov, linewidth=0, rasterized=True, cmap=cmap_cividis,
                                               norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=1., vmax=1000))
    divider = make_axes_locatable(ax_mean_markov)
    cax_mean_markov = divider.append_axes('right', size='5%', pad=0.05)
    cbar_mean_markov = fig.colorbar(cs_mean_markov, cax=cax_mean_markov, orientation='vertical')
    cbar_mean_markov.set_label(r"$\langle T \rangle$", loc="center")
    ax_mean_markov.set_title("Two-component")
    ax_mean_markov.set_ylabel(r"$p$")
    ax_mean_markov.set_xscale("log")
    ax_mean_markov.set_yscale("log")
    ax_mean_markov.set_xticklabels([])
    ax_cv_markov.set_ylabel(r"$p$")
    ax_cv_markov.set_xlabel(r"$\tau$")
    ax_cv_markov.set_xscale("log")
    ax_cv_markov.set_yscale("log")


    cs_cv_markov = ax_cv_markov.pcolormesh(taus, js, CV_markov, linewidth=0, rasterized=True, vmin=0., vmax=1.0, cmap=cmap_cividis)
    divider = make_axes_locatable(ax_cv_markov)
    cax_cv_markov = divider.append_axes('right', size='5%', pad=0.05)
    cbar_cv_markov = fig.colorbar(cs_cv_markov, cax=cax_cv_markov, orientation='vertical')
    cbar_cv_markov.set_label(r"$CV_T$", loc="center")
    cbar_cv_markov.set_ticks([0, 0.5, 1])


    for i, interpretation in enumerate(["Ito", "Stratonovich", "Hänggi"]):
        ax_mean = fig.add_axes([w0 + w + w_cb + dw1 + i*(w + dw2), h0 + h + dh, w, h])
        ax_cv = fig.add_axes([w0 + w + w_cb + dw1 + i*(w + dw2), h0, w, h])
        ax_mean.set_xscale("log")
        ax_mean.set_yscale("log")

        ax_cv.set_xlabel(r"$\tau$")
        ax_cv.set_xscale("log")
        ax_cv.set_yscale("log")
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
        elif interpretation == "Hänggi":
            interpret = "haenggi"
            ax_mean.set_xticklabels([])
            ax_mean.set_yticklabels([])
            ax_cv.set_yticklabels([])
            a = 1.0


        data_langevin = np.loadtxt(home + f"/Data/calcium_spikes_theory/langevin_{interpret:s}_ca_mean_CV_K{N:d}_N{n:d}_no_adap.dat")

        taus_m, js_m, Ts_m, cvs_m, num_m = np.transpose(data_markov)
        taus_l, js_l, Ts_l, cvs_l, num_l = np.transpose(data_langevin)

        taus = np.logspace(-1, 2, size)
        js = np.logspace(-3, 0, size)


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


        cs = ax_mean.pcolormesh(taus, js, dif_ISI, linewidth=0, rasterized=True, vmin=-0.4, vmax=0.4,
                                      cmap=cmap_coolwarm)
        # , norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))

        if interpretation == "Hänggi":
            divider = make_axes_locatable(ax_mean)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cbar = fig.colorbar(cs, cax=cax, orientation='vertical')
            cbar.set_ticks([-0.4, 0., 0.4])
            cbar.set_label(r"$\delta \langle T \rangle$", loc="center")

        ax_mean.set_title(f"{interpretation:s}")

        # ax2.set_ylim([0.01, 1])
        # ax2.set_xlim([0.1, 10])


        cs = ax_cv.pcolormesh(taus, js, dif_CV, linewidth=0, rasterized=True, vmin=-0.4, vmax=0.4, cmap=cmap_coolwarm)
        # , norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))

        if interpretation == "Hänggi":
            divider = make_axes_locatable(ax_cv)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cbar = fig.colorbar(cs, cax=cax, orientation='vertical')
            cbar.set_ticks([-0.4, 0., 0.4])
            cbar.set_label(r"$\delta CV_T$", loc="center")

        #ax4.set_ylim([0.01, 1])
        #ax4.set_xlim([0.1, 10])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig5.pdf", transparent=True)
    plt.show()


