import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import styles as st

if __name__ == "__main__":
    taumax = 1000
    taumin = 10
    epsmax = 1.00
    epsmin = 0.01
    tau_ers = np.logspace(1, 3, 41)
    eps_ers = np.logspace(-2, 0, 41)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(7, 5))
    gs = gridspec.GridSpec(2, 2)
    ax_cv1 = fig.add_subplot(gs[0, 0])
    ax_dcv1 = fig.add_subplot(gs[0, 1])
    ax_cv2 = fig.add_subplot(gs[1, 0])
    ax_dcv2 = fig.add_subplot(gs[1, 1])
    ax_cv1.text(0.1, 0.1, r"A$_1$", fontsize=13, transform=ax_cv1.transAxes)
    ax_dcv1.text(0.1, 0.1, r"A$_2$", fontsize=13, transform=ax_dcv1.transAxes)
    ax_cv2.text(0.1, 0.1, r"B$_1$", fontsize=13, c="w", transform=ax_cv2.transAxes)
    ax_dcv2.text(0.1, 0.1, r"B$_2$", fontsize=13, transform=ax_dcv2.transAxes)
    axis_cv = [ax_cv1, ax_cv2]
    axis_dcv = [ax_dcv1, ax_dcv2]
    ax_cv1.set_title("mean-driven")
    ax_cv2.set_title("excitable")

    ax_dcv1.set_title("mean-driven")
    ax_dcv2.set_title("excitable")

    ax_dcv1.scatter(500, 0.04, marker="x", c="k", zorder=3)
    ax_dcv1.scatter(500, 0.4, marker="x", c="k", zorder=3)
    ax_dcv2.scatter(500, 0.1, marker="x", c="k", zorder=3)
    colors = st.Colors()

    size = 41
    N = 10
    n = 5
    m = 3
    home = os.path.expanduser("~")
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    for ax_cv, ax_dcv, tau, p in zip(axis_cv, axis_dcv, taus, ps):
        ISIs_no_adap = np.loadtxt(home + f"/Data/calcium_spikes_markov/Data_no_adap/spike_times_markov_ip1.00_tau{tau:.2e}_j{p:.2e}_K10_5.dat")
        T0 = np.mean(ISIs_no_adap)
        CV0 = np.std(ISIs_no_adap)/T0

        data_markov = np.loadtxt(home + f"/Data/calcium_spikes_theory/markov_ca_stationary_statistics_tau{tau:.2e}_j{p:.2e}_K10_N5_adap.dat")
        taus_a, eps_a, Ts, cvs, count = np.transpose(data_markov)
        taus_a = taus_a.reshape(size, size)
        eps_a = eps_a.reshape(size, size)
        Ts = Ts.reshape(size, size)
        cvs = cvs.reshape(size, size)
        dcvs = (cvs - CV0)/CV0
        #ISIs = np.empty([size, size])
        #CVs = np.empty([size, size])
        #dCVs = np.empty([size, size])
        #for k, T in enumerate(Ts):
        #    ISIs[k // size, k % size] = T
        #for k, cv in enumerate(cvs):
        #    CVs[k // size, k % size] = cv
        #    dCVs[k // size, k % size] = (cv - CV0) / CV0

        cmap_cividis = plt.get_cmap("YlGnBu", 11)
        cmap_coolwarm = plt.get_cmap("coolwarm", 11)
        cmap_coolwarm_list = [cmap_coolwarm(i) for i in range(cmap_coolwarm.N)]
        cmap_coolwarm_list[5] = (1., 1., 1., 1.)
        cmap_coolwarm = mpl.colors.LinearSegmentedColormap.from_list('custom cmap', cmap_coolwarm_list, cmap_coolwarm.N)

        ax_cv.set_ylabel(r"$\varepsilon$")
        ax_cv.set_xlabel(r"$\tau_{er}$")
        ax_cv.set_xscale("log")
        ax_cv.set_yscale("log")
        ax_cv.set_ylim([epsmin, epsmax])
        ax_cv.set_xlim([taumin, taumax])
        ax_dcv.set_ylabel(r"$\varepsilon$")
        ax_dcv.set_xlabel(r"$\tau_{er}$")
        ax_dcv.set_xscale("log")
        ax_dcv.set_yscale("log")
        ax_dcv.set_ylim([epsmin, epsmax])
        ax_dcv.set_xlim([taumin, taumax])

        cs_cv_markov = ax_cv.pcolormesh(taus_a, eps_a, cvs, linewidth=0, rasterized=True, shading='gouraud', vmin=0.0, vmax=1.0, cmap=cmap_cividis)
        divider = make_axes_locatable(ax_cv)
        cax_cv_markov = divider.append_axes('right', size='5%', pad=0.05)
        cbar_cv_markov = fig.colorbar(cs_cv_markov, cax=cax_cv_markov, orientation='vertical')
        cbar_cv_markov.set_label(r"$CV_T$", loc="center")
        cbar_cv_markov.set_ticks([0, 0.5, 1.0])

        cs_dcv_markov = ax_dcv.pcolormesh(taus_a, eps_a, dcvs, linewidth=0, rasterized=True, shading='gouraud', vmin=-1.0, vmax=1.0, cmap=plt.get_cmap("RdBu_r", 11))
        #ax_cv.contour(taus_adap, epss_adap, CVs, levels=[0])
        divider = make_axes_locatable(ax_dcv)
        cax_dcv_markov = divider.append_axes('right', size='5%', pad=0.05)
        cbar_dcv_markov = fig.colorbar(cs_dcv_markov, cax=cax_dcv_markov, orientation='vertical')
        cbar_dcv_markov.set_label(r"$(CV_T - CV_{T,0})/CV_{T,0}$", loc="center")
        cbar_dcv_markov.set_ticks([-1., 0, 1.0])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig6.pdf", transparent=True)
    plt.show()


