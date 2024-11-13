import numpy as np
import matplotlib.pyplot as plt
import scipy
from  matplotlib import gridspec
import matplotlib.ticker as ticker
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    w = 3.25 * 1.25
    h = 2.50 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(nrows=2, ncols=2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])

    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"C", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"D", fontsize=11, transform=ax4.transAxes, va='top')

    home = os.path.expanduser("~")

    ax1.set_ylabel(r"$\langle T \rangle$")
    ax1.set_xlabel(r"$\lambda_{\rm ref}$")
    ax1.set_xlim([10, 100])
    ax1.set_ylim([0, 130])
    ax1.set_xscale("log")
    ax1.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax1.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax2.set_ylabel(r"$CV_T$")
    ax2.set_xlabel(r"$\lambda_{\rm ref}$")
    ax2.set_xscale("log")
    ax2.set_xlim([10, 100])
    ax2.set_ylim([0, 0.65])
    ax2.set_yticks([0, 0.25, 0.5])
    ax2.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax2.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax3.set_ylabel(r"$\mu_{\rm puff}$", fontsize=11)
    ax3.set_xticks([0.2, 0.5])
    ax3.set_ylim([0, 1.])
    ax3.set_xticklabels(["$c_0$", "$c_T$"])

    ax4.set_ylabel(r"$D_{\rm puff}$", fontsize=11)
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
        r_opn_single = 0.1
        mus = []
        Ds = []
        for ca in cas:
            mu = fc.mean_jp_single_theory(ca, N, M, 1, r_ref)
            D = fc.noise_intensity_jp_single_theory(ca, N, M, 1, r_ref = r_ref)
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
    ax3.legend(fancybox=False, fontsize=7, loc=1)


    home = os.path.expanduser("~")
    r_refs = np.logspace(1, 2, 11)
    tau = 5.
    j = 0.015
    mean_isis = []
    cv_isis = []
    mean_isis_langevin = []
    cv_isis_langevin = []
    r_opn_max = 0.
    for i, r_ref in enumerate(r_refs):
        folder = home + "/Data/calcium/markov/renewal/over_rref/"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_{i:d}.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        mean = np.mean(data_isis)
        cv = np.std(data_isis)/mean
        mean_isis.append(mean)
        cv_isis.append(cv)

        folder_langevin = home + "/Data/calcium/langevin/renewal/over_rref/"
        file_spikes_langevin = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_{i:d}.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_langevin = np.mean(data_isis_langevin)
        cv_langevin = np.std(data_isis_langevin)/mean_langevin
        mean_isis_langevin.append(mean_langevin)
        cv_isis_langevin.append(cv_langevin)


    ax1.scatter(np.append(r_refs[1:3], r_refs[4:]), mean_isis[1:3] + mean_isis[4:], s=20, fc="w", ec=st.colors[0], zorder=1)
    ax2.scatter(np.append(r_refs[1:3], r_refs[4:]), cv_isis[1:3] + cv_isis[4:], s=20, fc="w", ec=st.colors[0], zorder=1, label="Two-comp.")
    ax1.plot(r_refs[1:], mean_isis_langevin[1:], c=st.colors[0], zorder=0)
    ax2.plot(r_refs[1:], cv_isis_langevin[1:], c=st.colors[0], zorder=0, label="Langevin")
    ax1.scatter(r_refs[3], mean_isis[3], marker="X", s=40, fc="w", ec=st.colors[0])
    ax2.scatter(r_refs[3], cv_isis[3], marker="X", s=40, fc="w", ec=st.colors[0])
    ax2.legend(fancybox=False, fontsize=7, loc=1)

    def tau_p_bif(x):
        r_cls = 50
        r_opn_single = 0.1
        K = 10
        N = 5
        M = 3
        tau = 5
        p = 0.015

        ropnmax = N * r_opn_single * ((1. + np.power(0.20, 3)) / np.power(0.20, 3))
        r_opn_ct = ropnmax * (np.power(0.5, 3) / (1. + np.power(0.5, 3)))

        mean_puff = (6) * (7) / (6 * r_cls)
        tau_tot = 1 / r_opn_ct + (M-1) / x + (N+1) / (2 * r_cls)
        mean_puff_ct = mean_puff / tau_tot
        return (0.5 - 0.2) / (K * mean_puff_ct * tau) - p

    r_ref_bif = scipy.optimize.fsolve(tau_p_bif, x0 = np.array(20))

    ax1.axvspan(10, r_ref_bif, alpha=0.3, color="C7")
    ax1.axvline(r_ref_bif, lw=1, color="C7", zorder=1)

    ax2.axvspan(10, r_ref_bif, alpha=0.3, color="C7")
    ax2.axvline(r_ref_bif, lw=1, color="C7", zorder=1)
    #plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/SUB2/figures/fig9.pdf", transparent=True)
    plt.show()