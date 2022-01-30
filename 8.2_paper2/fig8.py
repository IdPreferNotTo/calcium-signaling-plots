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
    data_markov_adap = np.loadtxt(home + "/Data/Calcium/data/markov_ca_mean_CV_N{:d}_n{:d}_fix_adap_taua5.00e+02_ampa5.00e-02.dat".format(10, 5))

    taus_n, js_n, Ts_n, cvs_n, num_n = np.transpose(data_markov)
    taus_a, js_a, Ts_a, cvs_a, num_a = np.transpose(data_markov_adap)

    taus_n = np.logspace(-1, 2, 50)
    js_n = np.logspace(-3, 0, 50)

    taus_a = np.logspace(-1, 2, 31)
    js_a = np.logspace(-3, 0, 31)

    ISI_no_adap = np.empty([50, 50])
    CV_no_adap = np.empty([50, 50])
    ISI_adap = np.empty([31, 31])
    CV_adap = np.empty([31, 31])


    for n, T in enumerate(Ts_n):
        if T >= 1_000:
            ISI_no_adap[n // 50, n % 50] = np.nan
        else:
            ISI_no_adap[n // 50, n % 50] = T
    for n, cv in enumerate(cvs_n):
        if cv == 1:
            CV_no_adap[n // 50, n % 50] = np.nan
        else:
            CV_no_adap[n // 50, n % 50] = cv

    for n, T in enumerate(Ts_a):
        if T >= 1_000:
            ISI_adap[n // 31, n % 31] = np.nan
        else:
            ISI_adap[n // 31, n % 31] = T
    for n, cv in enumerate(cvs_a):
        if cv == 1:
            CV_adap[n // 31, n % 31] = np.nan
        else:
            CV_adap[n // 31, n % 31] = cv

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
    cs1 = ax1.pcolormesh(taus_n, js_n, ISI_no_adap, linewidth=0, rasterized=True, cmap=cmap_cividis, norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=1., vmax=1000))
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='5%', pad=0.05)
    cbar1 = fig.colorbar(cs1, cax=cax1, orientation='vertical')
    ax1.set_title(r"$\langle T \rangle$")
    #cbar1.set_ticks([0, 100, 200])
    ax1.set_ylabel(r"$P^*$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    #ax1.set_ylim([0.01, 1])
    #ax1.set_xlim([0.1, 10])

    cs2 = ax2.pcolormesh(taus_a, js_a, ISI_adap, linewidth=0, rasterized=True, cmap=cmap_cividis, norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=1., vmax=1000))
    # , norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes('right', size='5%', pad=0.05)
    cbar2 = fig.colorbar(cs2, cax=cax2, orientation='vertical')
    ax2.set_title(r"$\langle T \rangle$")
    #cbar1.set_ticks([0, 100, 200])
    ax2.set_ylabel(r"$P^*$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    cs3 = ax3.pcolormesh(taus_n, js_n, CV_no_adap, linewidth=0, rasterized=True, vmin=0., vmax=1.0, cmap=cmap_cividis)

    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes('right', size='5%', pad=0.05)
    cbar3 = fig.colorbar(cs3, cax=cax3, orientation='vertical')
    ax3.set_title('$C_V$')
    cbar3.set_ticks([0, 0.5, 1])
    ax3.set_ylabel(r"$P^*$")
    ax3.set_xlabel(r"$\tau$")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    #ax3.set_ylim([0.01, 1])
    #ax3.set_xlim([0.1, 100])

    cs4 = ax4.pcolormesh(taus_a, js_a, CV_adap, linewidth=0, rasterized=True, vmin=0., vmax=1.0, cmap=cmap_cividis)
    # , norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))

    divider = make_axes_locatable(ax4)
    cax4 = divider.append_axes('right', size='5%', pad=0.05)
    cbar4 = fig.colorbar(cs4, cax=cax4, orientation='vertical')
    ax4.set_title('$C_V$')
    cbar4.set_ticks([0, 0.5, 1])
    ax4.set_ylabel(r"$P^*$")
    ax4.set_xlabel(r"$\tau$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    #ax3.set_ylim([0.01, 1])
    #ax3.set_xlim([0.1, 100])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig8.png", transparent=True)
    plt.show()


