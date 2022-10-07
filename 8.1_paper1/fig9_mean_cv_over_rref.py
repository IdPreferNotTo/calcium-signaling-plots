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

    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    ax1.text(0.20, 0.95, r"A$_1$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.20, 0.95, r"B$_1$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.10, 0.95, r"A$_2$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.10, 0.95, r"B$_2$", fontsize=13, transform=ax4.transAxes, va='top')

    home = os.path.expanduser("~")

    ax1.set_ylabel(r"$\langle T \rangle$")
    ax1.set_xlabel(r"$\lambda_{\rm ref}$")
    ax1.set_xlim([10, 100])
    ax1.set_ylim([0, 100])
    ax1.set_xscale("log")
    ax1.axvspan(10, 20, alpha=0.3, color="C7")
    ax1.axvline(20, lw=1, color="C7", zorder=1)

    ax2.set_ylabel(r"$CV_T$")
    ax2.set_xlabel(r"$\lambda_{\rm ref}$")
    ax2.set_xscale("log")
    ax2.set_xlim([10, 100])
    ax2.set_ylim([0, 1])
    ax2.axvspan(10, 20, alpha=0.3, color="C7")
    ax2.axvline(20, lw=1, color="C7", zorder=1)

    ax3.set_ylabel(r"$\mu_{\rm puff}(c_{\rm i})$", fontsize=11)
    ax3.set_xticks([0.2, 0.5])
    ax3.set_ylim([0, 1.])
    ax3.set_xticklabels(["$c_0$", "$c_T$"])

    ax4.set_ylabel(r"$2D_{\rm puff}(c_{\rm i})$", fontsize=11)
    ax4.set_ylim([0, 0.07])
    ax4.set_xticks([0.2, 0.5])
    ax4.set_xticklabels(["$c_0$", "$c_T$"])

    N = 5
    M = 3
    mus = []
    Ds = []
    cas = np.linspace(0.2, 0.5, 100)
    r_refs = [1, 10, 100, 10**5]
    c0 = 0.2
    ropn_at_c0 = 0.5
    ropn_max = ropn_at_c0 * (1 + np.power(c0, 3))/np.power(c0, 3)
    print(ropn_max)
    for i, r_ref in enumerate(r_refs):
        r_opn_single = 1 / (N * (2.1 - (M - 1) / r_ref))
        mus = []
        Ds = []
        for ca in cas:
            mu = fc.mean_puff_single(ca, N, M, 1, r_opn_single, r_ref)
            D = fc.intensity_puff_single(ca, N, M, 1, r_opn_single, r_ref)
            mus.append(mu)
            Ds.append(D)
        if i == 3:
            ax3.plot(cas, mus, lw=1, color="k", zorder=1)
            ax4.plot(cas, Ds, lw=1, color="k", zorder=1)
        elif i == 2:
            ax3.plot(cas, mus, lw=1, color=st.colors[i + 2], label=rf"$\lambda_{{\rm ref}} = {r_ref:.0f}$")
            ax4.plot(cas, Ds, lw=1, color=st.colors[i + 2])
        else:
            ax3.plot(cas, mus, lw=1, color=st.colors[i + 1], label=rf"$\lambda_{{\rm ref}} = {r_ref:.0f}$")
            ax4.plot(cas, Ds, lw=1, color=st.colors[i + 1])
    ax3.legend(fancybox=False, fontsize=8, loc=9)


    home = os.path.expanduser("~")
    r_refs = np.logspace(1, 3, 21)
    taus = [5.62, 1.78]
    js = [0.0126, 0.0355]
    tau = taus[0]
    j = js[0]
    mean_isis = []
    cv_isis = []
    mean_isis_langevin = []
    cv_isis_langevin = []
    r_opn_max = 0.
    for i, r_ref in enumerate(r_refs):
        folder = home + "/Data/calcium_spikes_markov/Data_no_adap_rrefs/"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_{i:d}.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        mean = np.mean(data_isis)
        cv = np.std(data_isis)/mean
        mean_isis.append(mean)
        cv_isis.append(cv)

        folder_langevin = home + "/Data/calcium_spikes_langevin_strat/Data_no_adap_rrefs/"
        file_spikes_langevin = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_{i:d}.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_langevin = np.mean(data_isis_langevin)
        cv_langevin = np.std(data_isis_langevin)/mean_langevin
        mean_isis_langevin.append(mean_langevin)
        cv_isis_langevin.append(cv_langevin)


    ax1.plot(r_refs, mean_isis, color=st.colors[0], zorder=1, )
    ax2.plot(r_refs, cv_isis, color=st.colors[0], zorder=1, label="Two-component")
    ax1.scatter(r_refs, mean_isis_langevin, s=20, fc="w", ec=st.colors[0])
    ax2.scatter(r_refs, cv_isis_langevin, s=20, fc="w", ec=st.colors[0], label="Langevin")
    ax2.legend(fancybox=False, fontsize=8)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig9.pdf", transparent=True)
    plt.show()