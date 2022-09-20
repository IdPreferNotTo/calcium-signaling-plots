import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    fig = plt.figure(tight_layout=True, figsize=(6, 3.5))
    gs = gridspec.GridSpec(nrows=2, ncols=2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    # axin1 = ax1.inset_axes([0.45, 0.50, 0.50, 0.50])
    # axin2 = ax2.inset_axes([0.45, 0.50, 0.50, 0.50])
    # axin3 = ax3.inset_axes([0.55, 0.55, 0.45, 0.45])
    # axin4 = ax4.inset_axes([0.55, 0.55, 0.45, 0.45])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    ax1.text(0.05, 0.95, r"A$_{\rm i}$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.12, 0.95, r"A$_{\rm ii}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"B$_{\rm i}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.12, 0.95, r"B$_{\rm ii}$", fontsize=13, transform=ax4.transAxes, va='top')
    home = os.path.expanduser("~")


    ax1.set_ylabel(r"$\langle T \rangle$")
    ax1.set_xlabel(r"$K$")
    ax1.set_xlim([8.5, 15.5])
    ax1.set_ylim([0, 100])
    ax1.set_xticks([9, 10, 11, 12, 13, 14, 15])

    ax2.set_ylabel(r"$CV_T$")
    ax2.set_xlabel(r"$K$")
    ax2.set_xlim([8.5, 15.5])
    ax2.set_ylim([0, 0.5])
    ax2.set_xticks([9, 10, 11, 12, 13, 14, 15])

    ax3.set_ylabel(r"$\langle T \rangle$")
    ax3.set_xlabel(r"$K$")
    ax3.set_xlim([3.5, 10.5])
    ax3.set_ylim([0, 100])
    ax3.set_xticks([4, 5, 6, 7, 8, 9, 10])

    ax4.set_ylabel(r"$CV_T$")
    ax4.set_xlabel(r"$K$")
    ax4.set_xlim([3.5, 10.5])
    ax4.set_ylim([0, 0.5])
    ax4.set_xticks([4, 5, 6, 7, 8, 9, 10])

    # for axin in [axin1, axin3]:
    #     axin.set_title(r"$\mu_{\rm puff}(c_{\rm i})$", fontsize=11)
    #     axin.set_yticks([])
    #     axin.set_xticks([0.2, 0.5])
    #     axin.set_xticklabels(["$c_0$", "$c_T$"])
    #     axin.axhline(0, color="C7")

    # for axin in [axin2, axin4]:
    #     axin.set_title(r"$2D_{\rm puff}(c_{\rm i})$", fontsize=11)
    #     axin.set_yticks([])
    #     axin.set_xticks([0.2, 0.5])
    #     axin.set_xticklabels(["$c_0$", "$c_T$"])
    #     axin.text(0.95, 0.8, r"$K = 1$", color=st.colors[1], fontsize=8, transform=axin.transAxes,
    #                ha='right', va='top')
    #     axin.text(0.95, 0.55, r"$K = 5$", color=st.colors[2], fontsize=8, transform=axin.transAxes,
    #                ha='right', va='top')
    #     axin.text(0.95, 0.37, r"$K= 9$", color=st.colors[4], fontsize=8, transform=axin.transAxes,
    #                ha='right', va='top')

    N = 5
    M = 3
    tau = 5.62
    j = 0.0126
    r_ref = 20
    r_opn_single = 0.1

    mus = []
    Ds = []
    cas = np.linspace(0.2, 0.5, 100)
    Ks = [9, 12, 15]

    # for i, K in enumerate(Ks):
    #     mus = []
    #     Ds = []
    #     for ca in cas:
    #         mu = j*K*fc.mean_puff_single(ca, N, M, 1, r_opn_single, r_ref)
    #         D = j*j*K*fc.intensity_puff_single(ca, N, M, 1, r_opn_single, r_ref)
    #         mus.append(mu)
    #         Ds.append(D)
    #     if K == 9:
    #         axin1.plot(cas, mus, lw=1, color=st.colors[i + 2])
    #         axin2.plot(cas, Ds, lw=1, color=st.colors[i + 2])
    #     else:
    #         axin1.plot(cas, mus, lw=1, color=st.colors[i + 1])
    #         axin2.plot(cas, Ds, lw=1, color=st.colors[i + 1])

    mean_isis_langevin = []
    cv_isis_langevin = []
    mean_isis = []
    cv_isis = []
    Ks = np.arange(9, 16)
    for i, K in enumerate(Ks):
        folder = home + "/Data/calcium_spikes_markov/Data_no_adap/"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K{K:d}_5.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        print(len(data_isis))
        mean = np.mean(data_isis)
        cv = np.std(data_isis)/mean
        mean_isis.append(mean)
        cv_isis.append(cv)

        folder_langevin = home + "/Data/calcium_spikes_langevin_strat/Data_no_adap/"
        file_spikes_langevin = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{j:.2e}_K{K:d}_5.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_langevin = np.mean(data_isis_langevin)
        cv_langevin = np.std(data_isis_langevin)/mean_langevin
        mean_isis_langevin.append(mean_langevin)
        cv_isis_langevin.append(cv_langevin)

    ax1.plot(Ks, mean_isis, color=st.colors[0], zorder=1)
    ax2.plot(Ks, cv_isis, color=st.colors[0], zorder=1, label="Two-component")
    ax1.scatter(Ks, mean_isis_langevin, s=20, fc="w", ec=st.colors[0])
    ax1.text(0.95, 0.95, r"$p = {\rm const.}$", fontsize=11, transform=ax1.transAxes, ha="right", va='top')
    ax2.scatter(Ks, cv_isis_langevin, s=20, fc="w", ec=st.colors[0], label="Langevin")
    ax2.legend(fancybox=False, fontsize=8)


    mean_isis = []
    cv_isis = []
    mean_isis_langevin = []
    cv_isis_langevin = []

    N = 5
    M = 3
    mus = []
    Ds = []
    cas = np.linspace(0.2, 0.5, 100)
    Ks = [1, 5, 9]
    j = 0.126
    r_ref = 20
    c0 = 0.2
    # for i, K in enumerate(Ks):
    #     jsingle = j / K
    #     mus = []
    #     Ds = []
    #     for ca in cas:
    #         mu = jsingle * K * fc.mean_puff_single(ca, N, M, 1, r_opn_single, r_ref)
    #         D = jsingle * jsingle * K * fc.intensity_puff_single(ca, N, M, 1, r_opn_single, r_ref)
    #         mus.append(mu)
    #         Ds.append(D)
    #     if K == 9:
    #         axin3.plot(cas, mus, lw=1, color=st.colors[i + 2])
    #         axin4.plot(cas, Ds, lw=1, color=st.colors[i + 2])
    #     else:
    #         axin3.plot(cas, mus, lw=1, color=st.colors[i + 1])
    #         axin4.plot(cas, Ds, lw=1, color=st.colors[i + 1])

    mean_isis_langevin = []
    cv_isis_langevin = []
    mean_isis = []
    cv_isis = []
    Ks = np.arange(4, 11)
    for i, K in enumerate(Ks):
        jsingle = j / K
        folder = home + "/Data/calcium_spikes_markov/Data_no_adap/"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{jsingle:.2e}_K{K:d}_5.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        print(len(data_isis))
        mean = np.mean(data_isis)
        cv = np.std(data_isis) / mean
        mean_isis.append(mean)
        cv_isis.append(cv)

        folder_langevin = home + "/Data/calcium_spikes_langevin_strat/Data_no_adap/"
        file_spikes_langevin = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{jsingle:.2e}_K{K:d}_5.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_langevin = np.mean(data_isis_langevin)
        cv_langevin = np.std(data_isis_langevin) / mean_langevin
        mean_isis_langevin.append(mean_langevin)
        cv_isis_langevin.append(cv_langevin)

    ax3.plot(Ks, mean_isis, color=st.colors[0], zorder=1)
    ax4.plot(Ks, cv_isis, color=st.colors[0], zorder=1)
    ax3.scatter(Ks, mean_isis_langevin, s=20, fc="w", ec=st.colors[0])
    ax3.text(0.95, 0.95, r"$pK = {\rm const.}$", fontsize=11, transform=ax3.transAxes, ha="right", va='top')
    ax4.scatter(Ks, cv_isis_langevin, s=20, fc="w", ec=st.colors[0])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig10.pdf", transparent=True)
    plt.show()