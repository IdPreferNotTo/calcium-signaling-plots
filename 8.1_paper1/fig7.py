import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, cm
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import styles as st

if __name__ == "__main__":
    home = os.path.expanduser("~")
    data_markov = np.loadtxt(home + "/Data/Calcium/data/markov_ca_mean_CV_N{:d}_n{:d}_no_adap.dat".format(10, 5))
    data_langevin = np.loadtxt(home + "/Data/Calcium/data/langevin_ca_mean_CV_N{:d}_n{:d}_no_adap.dat".format(10, 5))

    taus_m, js_m, Ts_m, cvs_m, num_m = np.transpose(data_markov)
    taus_l, js_l, Ts_l, cvs_l, num_l = np.transpose(data_langevin)

    taus = np.logspace(-1, 2, 50)
    js = np.logspace(-3, 0, 50)

    ISI_markov = np.empty([50, 50])
    CV_markov = np.empty([50, 50])
    ISI_langevin = np.empty([50, 50])
    CV_langevin = np.empty([50, 50])
    dif_ISI = np.empty([50, 50])
    dif_CV = np.empty([50, 50])

    for n, T in enumerate(Ts_l):
        if T >= 1_000:
            ISI_langevin[n // 50, n % 50] = np.nan
        else:
            ISI_langevin[n // 50, n % 50] = T
    for n, cv in enumerate(cvs_l):
        if cv == 1:
            CV_langevin[n // 50, n % 50] = np.nan
        else:
            CV_langevin[n // 50, n % 50] = cv

    for n, T in enumerate(Ts_m):
        if T >= 1_000:
            ISI_markov[n // 50, n % 50] = np.nan
            dif_ISI[n // 50, n % 50] = np.nan
        else:
            ISI_markov[n // 50, n % 50] = T
            dif_ISI[n // 50, n % 50] = -(T - ISI_langevin[n // 50, n % 50]) / T
    for n, cv in enumerate(cvs_m):
        if cv == 1:
            CV_markov[n // 50, n % 50] = np.nan
            dif_CV[n // 50, n % 50] = np.nan
        else:
            CV_markov[n // 50, n % 50] = cv
            dif_CV[n // 50, n % 50] = -(cv - CV_langevin[n // 50, n % 50]) / cv




    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 4.))
    grids = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(grids[0, 0])
    ax2 = fig.add_subplot(grids[0, 1])
    ax3 = fig.add_subplot(grids[1, 0])
    ax4 = fig.add_subplot(grids[1, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    cmap_cividis = plt.get_cmap("cividis", 8)
    cmap_coolwarm = plt.get_cmap("coolwarm", 8)
    cs1 = ax1.pcolormesh(taus, js, ISI_markov, linewidth=0, rasterized=True, cmap=cmap_cividis, norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=1., vmax=1000))
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='5%', pad=0.05)
    cbar1 = fig.colorbar(cs1, cax=cax1, orientation='vertical')
    ax1.set_title(r"$\langle T \rangle$")
    #cbar1.set_ticks([0, 100, 200])
    ax1.set_ylabel(r"$D^*_{\rm puff}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    #ax1.set_ylim([0.01, 1])
    #ax1.set_xlim([0.1, 10])

    cs2 = ax2.pcolormesh(taus, js, dif_ISI, linewidth=0, rasterized=True, vmin=-0.2, vmax=0.2, cmap=cmap_coolwarm)
    # , norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes('right', size='5%', pad=0.05)
    cbar2 = fig.colorbar(cs2, cax=cax2, orientation='vertical')
    ax2.set_title(r"$\delta \langle T \rangle$")
    cbar2.set_ticks([-0.2, 0., 0.2])
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    #ax2.set_ylim([0.01, 1])
    #ax2.set_xlim([0.1, 10])

    cs3 = ax3.pcolormesh(taus, js, CV_markov, linewidth=0, rasterized=True, vmin=0., vmax=1.0, cmap=cmap_cividis)

    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes('right', size='5%', pad=0.05)
    cbar3 = fig.colorbar(cs3, cax=cax3, orientation='vertical')
    ax3.set_title('$C_V$')
    cbar3.set_ticks([0, 0.5, 1])
    ax3.set_ylabel(r"$D^*_{\rm puff}$")
    ax3.set_xlabel(r"$\tau$")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    #ax3.set_ylim([0.01, 1])
    #ax3.set_xlim([0.1, 100])

    cs4 = ax4.pcolormesh(taus, js, dif_CV, linewidth=0, rasterized=True, vmin=-0.4, vmax=0.4, cmap=cmap_coolwarm)
    # , norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))

    divider = make_axes_locatable(ax4)
    cax4 = divider.append_axes('right', size='5%', pad=0.05)
    cbar4 = fig.colorbar(cs4, cax=cax4, orientation='vertical')
    ax4.set_title('$\delta C_V$')
    cbar4.set_ticks([-0.4, 0., 0.4])
    ax4.set_xlabel(r"$\tau$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    #ax4.set_ylim([0.01, 1])
    #ax4.set_xlim([0.1, 10])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig7.pdf", transparent=True)
    plt.show()


