import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import styles as st
import default_parameters as df
import functions as fc

if __name__ == "__main__":
    taumax = 1000
    taumin = 10
    epsmax = 1.00
    epsmin = 0.01
    tau_ers = np.logspace(1, 3, 41)
    eps_ers = np.logspace(-2, 0, 41)

    st.set_default_plot_style()
    w = 6.5 * 1.25
    h = 4 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(3, 3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])
    ax1.text(0.1, 0.1, r"A$_1$", fontsize=11, transform=ax1.transAxes)
    ax2.text(0.1, 0.1, r"A$_2$", fontsize=11, transform=ax2.transAxes)
    ax3.text(0.1, 0.1, r"A$_3$", fontsize=11, transform=ax3.transAxes)
    ax4.text(0.1, 0.1, r"B$_1$", fontsize=11, transform=ax4.transAxes)
    ax5.text(0.1, 0.1, r"B$_2$", fontsize=11, c="w", transform=ax5.transAxes)
    ax6.text(0.1, 0.1, r"B$_3$", fontsize=11, transform=ax6.transAxes)
    axis_cv = [ax2, ax5]
    axis_dcv = [ax3, ax6]
    gs1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs[2, 0], hspace=0.1)
    gs2 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs[2, 1], hspace=0.1)
    gs3 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs[2, 2], hspace=0.1)
    ax7 = fig.add_subplot(gs1[0])
    ax8 = fig.add_subplot(gs1[1])
    ax9 = fig.add_subplot(gs2[0])
    ax10 = fig.add_subplot(gs2[1])
    ax11 = fig.add_subplot(gs3[0])
    ax12 = fig.add_subplot(gs3[1])
    st.remove_top_right_axis([ax7, ax8, ax9, ax10, ax11, ax12])

    ax1.set_title("mean-driven")
    ax2.set_title("mean-driven")
    ax3.set_title("mean-driven")

    ax4.set_title("excitable")
    ax5.set_title("excitable")
    ax6.set_title("excitable")

    ax1.scatter(500, 0.02, marker="x", c="k", zorder=3)
    ax1.scatter(500, 0.2, marker="*", c="k", zorder=3)
    ax4.scatter(500, 0.1, marker="p", c="k", zorder=3)

    ax2.scatter(500, 0.02, marker="x", c="k", zorder=3)
    ax2.scatter(500, 0.2, marker="*", c="k", zorder=3)
    ax5.scatter(500, 0.1, marker="p", c="k", zorder=3)

    ax3.scatter(500, 0.02, marker="x", c="k", zorder=3)
    ax3.scatter(500, 0.2, marker="*", c="k", zorder=3)
    ax6.scatter(500, 0.1, marker="p", c="k", zorder=3)
    colors = st.Colors()

    size = 41
    N = 10
    n = 5
    m = 3
    home = os.path.expanduser("~")
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    for ax1, ax2, ax3, tau, p in zip([ax1, ax4], [ax2, ax5], [ax3, ax6], taus, ps):
        ax1.set_ylabel(r"$\varepsilon$")
        ax1.set_xlabel(r"$\tau_{\rm er}$")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylim([epsmin, epsmax])
        ax1.set_xlim([taumin, taumax])
        ax2.set_ylabel(r"$\varepsilon$")
        ax2.set_xlabel(r"$\tau_{\rm er}$")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_ylim([epsmin, epsmax])
        ax2.set_xlim([taumin, taumax])
        ax3.set_ylabel(r"$\varepsilon$")
        ax3.set_xlabel(r"$\tau_{\rm er}$")
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_ylim([epsmin, epsmax])
        ax3.set_xlim([taumin, taumax])

        ISIs_no_adap = np.loadtxt(home + f"/Data/calcium_spikes_markov/Data_no_adap/spike_times_markov_ip1.00_tau{tau:.2e}_j{p:.2e}_K10_5.dat")
        T0 = np.mean(ISIs_no_adap)
        CV0 = np.std(ISIs_no_adap)/T0

        data_markov = np.loadtxt(home + f"/Data/calcium/theory/markov_mean_cv_tau{tau:.2e}_j{p:.2e}_K10_N5_adap.dat")
        taus_a, eps_a, Ts, cvs, count = np.transpose(data_markov)
        taus_a = taus_a.reshape(size, size)
        eps_a = eps_a.reshape(size, size)
        Ts = Ts.reshape(size, size)
        cvs = cvs.reshape(size, size)
        dcvs = (cvs - CV0)/CV0


        cmap_cividis = plt.get_cmap("YlGnBu", 11)
        cmap_coolwarm = plt.get_cmap("coolwarm", 11)
        cmap_coolwarm_list = [cmap_coolwarm(i) for i in range(cmap_coolwarm.N)]
        cmap_coolwarm_list[5] = (1., 1., 1., 1.)
        cmap_coolwarm = mpl.colors.LinearSegmentedColormap.from_list('custom cmap', cmap_coolwarm_list, cmap_coolwarm.N)

        cs_T_markov = ax1.pcolormesh(taus_a, eps_a, Ts, linewidth=0, rasterized=True, shading='gouraud', norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=10., vmax=10_000), cmap=cmap_cividis)
        divider = make_axes_locatable(ax1)
        cax_T_markov = divider.append_axes('right', size='5%', pad=0.05)
        cbar_T_markov = fig.colorbar(cs_T_markov, cax=cax_T_markov, orientation='vertical')
        cbar_T_markov.set_label(r"$\langle T \rangle$", loc="center")
        #cbar_T_markov.set_ticks([0, 0.5, 1.0])

        cs_cv_markov = ax2.pcolormesh(taus_a, eps_a, cvs, linewidth=0, rasterized=True, shading='gouraud', vmin=0.0, vmax=1.0, cmap=cmap_cividis)
        divider = make_axes_locatable(ax2)
        cax_cv_markov = divider.append_axes('right', size='5%', pad=0.05)
        cbar_cv_markov = fig.colorbar(cs_cv_markov, cax=cax_cv_markov, orientation='vertical')
        cbar_cv_markov.set_label(r"$C_V$", loc="center")
        cbar_cv_markov.set_ticks([0, 0.5, 1.0])

        cs_dcv_markov = ax3.pcolormesh(taus_a, eps_a, dcvs, linewidth=0, rasterized=True, shading='gouraud', vmin=-1.0, vmax=1.0, cmap=plt.get_cmap("RdBu_r", 11))
        #ax_cv.contour(taus_adap, epss_adap, CVs, levels=[0])
        divider = make_axes_locatable(ax3)
        cax_dcv_markov = divider.append_axes('right', size='5%', pad=0.05)
        cbar_dcv_markov = fig.colorbar(cs_dcv_markov, cax=cax_dcv_markov, orientation='vertical')
        cbar_dcv_markov.set_label(r"$\delta C_V$", loc="center")
        cbar_dcv_markov.set_ticks([-1., 0, 1.0])

    ax7.text(0.1, 0.9, r"C", fontsize=11, transform=ax7.transAxes)
    ax9.text(0.1, 0.9, r"D", fontsize=11, transform=ax9.transAxes)
    ax11.text(0.1, 0.9, r"E", fontsize=11, transform=ax11.transAxes)
    ax7.scatter(2.5, 4., s=50, marker="*", c="k", zorder=3)
    ax9.scatter(2.5, 1.6, s=50, marker="x", c="k", zorder=3)
    ax11.scatter(2.5, 1.3, s=50, marker="p", c="k", zorder=3)
    #ax8.set_ylim([-5, 5])
    #ax10.set_ylim([-1.1, 1.1])
    #ax12.set_ylim([-1.1, 1.1])
    for ax1, ax2, tau, p, tau_er, eps_er in zip([ax7, ax9, ax11], [ax8, ax10, ax12], [5, 5, 1], [0.015, 0.015, 0.06], [500, 500, 500], [0.2, 0.02, 0.1]):
        ax1.set_xlim([0, 3])
        ax2.set_xlim([0, 3])
        ax1.set_xlabel("")
        ax2.set_xlabel(r"$t / \langle T \rangle$")
        ax1.set_xticks([1, 2, 3])
        ax2.set_xticks([1, 2, 3])
        ax1.set_xticklabels([])
        ax1.set_ylabel(r"$p_{\rm ISI}(t)$")
        ax2.set_ylabel(r"$C_z(t)$")

        data_isi = df.load_spike_times_markov(tau, p, cer=False)
        mean = np.mean(data_isi)
        s_times = fc.isis_to_spike_times(data_isi)
        ts, corr = fc.autocorrelation_from_spike_times(s_times, mean)
        ax1.hist(data_isi/mean, density=True, fc="C7", bins=25, alpha=0.75)
        ax2.plot(ts/mean, corr/abs(corr[0]), c="C7")

        data_isi = np.loadtxt(home + f"/Data/calcium/markov/adaptive/sequence/"
                                     f"long_spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat")
        print(data_isi[0])
        print(data_isi.size)
        mean = np.mean(data_isi)
        s_times = fc.isis_to_spike_times(data_isi)
        ts, corr = fc.autocorrelation_from_spike_times(s_times, mean)
        ax1.hist(data_isi/mean, density=True, fc=st.colors[1], bins=25, alpha=0.75)
        ax2.plot(ts/mean, corr/abs(corr[0]))


    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig4.pdf", transparent=True)
    plt.show()


