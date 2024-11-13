import numpy as np
import matplotlib.pyplot as plt
import scipy
from  matplotlib import gridspec
import os

import styles as st

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    w = 3.25 * 1.25
    h = 2.5 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(nrows=2, ncols=2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    ax1.text(0.05, 0.95, r"A$_1$", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"A$_2$", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"B$_1$", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"B$_2$", fontsize=11, transform=ax4.transAxes, va='top')
    home = os.path.expanduser("~")

    ax1.set_ylabel(r"$\langle T \rangle$")
    ax1.set_xlabel(r"$K$")
    ax1.set_xlim([7.5, 14.5])
    ax1.set_ylim([0, 100])
    ax1.set_xticks([8, 9, 10, 11, 12, 13, 14])

    ax2.set_ylabel(r"$CV_T$")
    ax2.set_xlabel(r"$K$")
    ax2.set_xlim([7.5, 14.5])
    ax2.set_ylim([0, 0.5])
    ax2.set_xticks([8, 9, 10, 11, 12, 13, 14])

    ax3.set_ylabel(r"$\langle T \rangle$")
    ax3.set_xlabel(r"$K$")
    ax3.set_xlim([7.5, 14.5])
    ax3.set_ylim([0, 100])
    ax3.set_xticks([8, 9, 10, 11, 12, 13, 14])

    ax4.set_ylabel(r"$CV_T$")
    ax4.set_xlabel(r"$K$")
    ax4.set_xlim([7.5, 14.5])
    ax4.set_ylim([0, 0.5])
    ax4.set_xticks([8, 9, 10, 11, 12, 13, 14])

    # MEAN-DRIVEN
    tau = 5.
    j = 0.015
    mean_isis_langevin = []
    cv_isis_langevin = []
    mean_isis = []
    cv_isis = []
    Ks = np.arange(9, 15)
    for i, K in enumerate(Ks):
        folder = home + "/Data/calcium/markov/renewal/over_K/"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K{K:d}_0.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        print(len(data_isis))
        mean = np.mean(data_isis)
        cv = np.std(data_isis)/mean
        mean_isis.append(mean)
        cv_isis.append(cv)

        folder_langevin = home + "/Data/calcium/langevin/renewal/over_K/"
        file_spikes_langevin = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{j:.2e}_K{K:d}_0.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_langevin = np.mean(data_isis_langevin)
        cv_langevin = np.std(data_isis_langevin)/mean_langevin
        mean_isis_langevin.append(mean_langevin)
        cv_isis_langevin.append(cv_langevin)

    ax1.scatter(np.append(Ks[0:1], Ks[2:]), np.append(mean_isis[0:1], mean_isis[2:]), s=20, fc="w", ec=st.colors[0], zorder=1)
    ax2.scatter(np.append(Ks[0:1], Ks[2:]), np.append(cv_isis[0:1], cv_isis[2:]), s=20, fc="w", ec=st.colors[0], zorder=1, label="Two-comp.")
    ax1.scatter(Ks[1], mean_isis[1], marker="X", s=40, fc="w", ec=st.colors[0])
    ax2.scatter(Ks[1], cv_isis[1], marker="X", s=40, fc="w", ec=st.colors[0])

    ax1.plot(Ks, mean_isis_langevin, c=st.colors[0], zorder=0)
    ax2.plot(Ks, cv_isis_langevin, c=st.colors[0], zorder=0, label="Langevin")
    ax1.text(0.95, 0.95, r"$p = {\rm const.}$", fontsize=11, transform=ax1.transAxes, ha="right", va='top')
    ax2.legend(fancybox=False, fontsize=7)

    tau = 5
    j = 0.15
    mean_isis_langevin = []
    cv_isis_langevin = []
    mean_isis = []
    cv_isis = []
    Ks = np.arange(8, 15)
    for i, K in enumerate(Ks):
        jsingle = j / K
        folder = home + "/Data/calcium/markov/renewal/over_K_and_p/"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{jsingle:.2e}_K{K:d}_0.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        print(len(data_isis))
        mean = np.mean(data_isis)
        cv = np.std(data_isis) / mean
        mean_isis.append(mean)
        cv_isis.append(cv)

        folder_langevin = home + "/Data/calcium/langevin/renewal/over_K_and_p/"
        file_spikes_langevin = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{jsingle:.2e}_K{K:d}_0.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_langevin = np.mean(data_isis_langevin)
        cv_langevin = np.std(data_isis_langevin) / mean_langevin
        mean_isis_langevin.append(mean_langevin)
        cv_isis_langevin.append(cv_langevin)

    ax3.scatter(np.append(Ks[0:2], Ks[3:]), np.append(mean_isis[0:1], mean_isis[2:]), s=20, fc="w", ec=st.colors[0], zorder=1)
    ax4.scatter(np.append(Ks[0:2], Ks[3:]), np.append(cv_isis[0:1], cv_isis[2:]), s=20, fc="w", ec=st.colors[0], zorder=1)
    ax3.scatter(Ks[2], mean_isis[2], marker="X", s=40, fc="w", ec=st.colors[0])
    ax4.scatter(Ks[2], cv_isis[2], marker="X", s=40, fc="w", ec=st.colors[0])

    ax3.plot(Ks, mean_isis_langevin, c=st.colors[0], zorder=0)
    ax4.plot(Ks, cv_isis_langevin, c=st.colors[0], zorder=0)
    ax3.text(0.95, 0.95, r"$pK = {\rm const.}$", fontsize=11, transform=ax3.transAxes, ha="right", va='top')


    def tau_p_bif(x):
        r_cls = 50
        r_ref = 20
        r_opn_single = 0.1
        N = 5
        M = 3
        tau = 5
        p = 0.015

        ropnmax = N * r_opn_single * ((1. + np.power(0.20, 3)) / np.power(0.20, 3))
        r_opn_ct = ropnmax * (np.power(0.5, 3) / (1. + np.power(0.5, 3)))

        mean_puff = (6) * (7) / (6 * r_cls)
        tau_tot = 1 / r_opn_ct + (M-1) / r_ref + (N+1) / (2 * r_cls)
        mean_puff_ct = mean_puff / tau_tot
        return (0.5 - 0.2) / (x * mean_puff_ct * tau) - p

    K_bif = scipy.optimize.fsolve(tau_p_bif, x0 = np.array(20))
    ax1.axvspan(7, K_bif, alpha=0.3, color="C7")
    ax1.axvline(K_bif, lw=1, color="C7", zorder=1)
    ax2.axvspan(7, K_bif, alpha=0.3, color="C7")
    ax2.axvline(K_bif, lw=1, color="C7", zorder=1)

    #plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/SUB2/figures/fig10.pdf", transparent=True)
    plt.show()