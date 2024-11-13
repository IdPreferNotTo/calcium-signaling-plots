import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df


if __name__ == "__main__":
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25*1.25
    h = 4.0*1.25
    # fig = plt.figure(tight_layout=True, figsize=(w, h))
    # gs = gridspec.GridSpec(3, 2)
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax2 = fig.add_subplot(gs[1, 0])
    # ax3 = fig.add_subplot(gs[2, 0])
    # ax4 = fig.add_subplot(gs[0, 1])
    # ax5 = fig.add_subplot(gs[1, 1])
    # ax6 = fig.add_subplot(gs[2, 1])

    fig, axs = plt.subplots(nrows=3, ncols=2, layout="constrained", figsize=(w, h))
    ax1 = axs[0, 0]
    ax2 = axs[1, 0]
    ax3 = axs[2, 0]
    ax4 = axs[0, 1]
    ax5 = axs[1, 1]
    ax6 = axs[2, 1]

    axis_top = [ax1, ax2, ax3]
    axis_bot = [ax4, ax5, ax6]
    axis = axis_top + axis_bot
    st.remove_top_right_axis(axis)
    ax1.text(0.15, 0.95, r"A$_1$", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.15, 0.95, r"A$_2$", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.15, 0.95, r"A$_3$", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.15, 0.95, r"B$_1$", fontsize=11, transform=ax4.transAxes, va='top')
    ax5.text(0.15, 0.95, r"B$_2$", fontsize=11, transform=ax5.transAxes, va='top')
    ax6.text(0.15, 0.95, r"B$_3$", fontsize=11, transform=ax6.transAxes, va='top')

    ax1.set_ylabel(r"$\langle T \rangle$ / s")
    ax1.set_xlabel(r"$\tau_{er}$ / s")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylim([20, 2000])

    ax2.set_ylabel(r"$C_V$")
    ax2.set_xlabel(r"$\tau_{er}$ / s")
    ax2.set_xscale("log")
    ax2.set_ylim([0, 1])

    ax3.set_ylabel(r"$\rho_1$")
    ax3.set_xlabel(r"$\tau_{er}$ / s")
    ax3.set_xscale("log")
    ax3.set_ylim([-0.5, 0.2])
    ax3.axhline(0, ls=":", c="k")

    ax4.set_ylabel(r"$\langle T \rangle$ / s")
    ax4.set_xlabel(r"$\varepsilon$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_ylim([20, 2000])

    ax5.set_ylabel(r"$C_V$")
    ax5.set_xlabel(r"$\varepsilon$")
    ax5.set_xscale("log")
    ax5.set_ylim([0, 1])

    ax6.set_ylabel(r"$\rho_1$")
    ax6.set_xlabel(r"$\varepsilon$")
    ax6.set_xscale("log")
    ax6.set_ylim([-0.5, 0.2])
    ax6.axhline(0, ls=":", c="k")
    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.03
    tau_er_fix = 300
    tau_ers = np.logspace(1, 3, 11)
    eps_ers = np.logspace(-2, 0, 11)

    #colors = ["#007A45", "#20639B"]
    colors = [st.colors[1], st.colors[3]]
    for i, (tau, p) in enumerate(zip(taus, ps)):
        if i == 0:
            label = "mean-driven"
        else:
            label = "excitable"
        means = []
        cvs = []
        p1s = []
        means_theory = []
        cv_theory = []
        means_langevin = []
        cvs_langevin = []
        p1s_langevin = []
        p2s_langevin = []
        for ii, tau_er in enumerate(tau_ers):
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er_fix)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            cvs.append(cv)
            p1s.append(p1)

            data_isi_langevin = df.load_spike_times_langevin(tau, p, cer=True, taua=tau_er, ampa=eps_er_fix)
            mean = np.mean(data_isi_langevin)
            std = np.std(data_isi_langevin)
            var = np.var(data_isi_langevin)
            cv = std / mean
            p1 = fc.k_corr(data_isi_langevin, data_isi_langevin, 1) / var
            p2 = fc.k_corr(data_isi_langevin, data_isi_langevin, 2) / var
            means_langevin.append(mean)
            cvs_langevin.append(cv)
            p1s_langevin.append(p1)
            p2s_langevin.append(p2)

            mean_isi_theory = fc.self_consistent_T_infty(tau, p, tau_er, eps_er_fix)
            means_theory.append(mean_isi_theory)
            cv_isi_theory = fc.self_consistent_CV_infty(tau, p, tau_er, eps_er_fix)
            cv_theory.append(cv_isi_theory)

        ax1.plot(tau_ers, means_langevin, c=colors[i], zorder=1)
        ax2.plot(tau_ers, cvs_langevin, c=colors[i], zorder=1)
        ax3.plot(tau_ers, p1s_langevin, c=colors[i], zorder=1)

        ax1.scatter(tau_ers, means, fc="w", ec=colors[i], s=20, label=label, zorder=2)
        ax2.scatter(tau_ers, cvs, fc="w", ec=colors[i], s=20, zorder=2)
        ax3.scatter(tau_ers, p1s, fc="w", ec=colors[i], s=20, zorder=2)

        if i == 1:
            # ax1.plot(tau_ers, means_theory, c=st.colors[7], ls=(0, (3, 1)), zorder=3, label="theory")
            ax1.plot(tau_ers, means_theory, c=colors[i], ls=(0, (3, 1)), zorder=3, label="theory")
        else:
            ax1.plot(tau_ers, means_theory, c=colors[i], ls=(0, (3, 1)), zorder=3)
            ax2.plot(tau_ers, cv_theory, c=colors[i], ls=(0, (3, 1)), zorder=3)
            # ax1.plot(tau_ers, means_theory, c=st.colors[7], ls=(0, (3, 1)), zorder=3)
            # ax2.plot(tau_ers, cv_theory, c=st.colors[7], ls=(0, (3, 1)), zorder=3)


        means = []
        cvs = []
        p1s = []
        means_theory = []
        cv_theory = []
        means_langevin = []
        cvs_langevin = []
        p1s_langevin = []

        for ii, eps_er in enumerate(eps_ers):
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er_fix, ampa=eps_er)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            cvs.append(cv)
            p1s.append(p1)

            data_isi_langevin = df.load_spike_times_langevin(tau, p, cer=True, taua=tau_er_fix, ampa=eps_er)
            mean = np.mean(data_isi_langevin)
            std = np.std(data_isi_langevin)
            var = np.var(data_isi_langevin)
            cv = std / mean
            p1 = fc.k_corr(data_isi_langevin, data_isi_langevin, 1) / var
            means_langevin.append(mean)
            cvs_langevin.append(cv)
            p1s_langevin.append(p1)

            mean_isi_theory = fc.self_consistent_T_infty(tau, p, tau_er_fix, eps_er)
            means_theory.append(mean_isi_theory)

            cv_isi_theory = fc.self_consistent_CV_infty(tau, p, tau_er_fix, eps_er)
            cv_theory.append(cv_isi_theory)


        ax4.plot(eps_ers, means_langevin, c=colors[i], zorder=1)
        ax5.plot(eps_ers, cvs_langevin, c=colors[i], zorder=1)
        ax6.plot(eps_ers, p1s_langevin, c=colors[i], zorder=1)

        ax4.scatter(eps_ers, means, fc="w", ec=colors[i], s=20, zorder=2)
        ax5.scatter(eps_ers, cvs, fc="w", ec=colors[i], s=20, zorder=2)
        ax6.scatter(eps_ers, p1s, fc="w", ec=colors[i], s=20, zorder=2)

        ax4.plot(eps_ers, means_theory, c=colors[i], ls=(0, (3, 1)), zorder=3)
        # ax4.plot(eps_ers, means_theory, c=st.colors[7], ls=(0, (3, 1)), zorder=3)
        if i == 0:
            ax5.plot(eps_ers, cv_theory, c=colors[i], ls=(0, (3, 1)), zorder=3)
            # ax5.plot(eps_ers, cv_theory, c=st.colors[7], ls=(0, (3, 1)), zorder=3)

    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.45, .0), loc=3,
                      ncol=4, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(.75)
    leg.legendHandles[2].set_color("k")
    home = os.path.expanduser("~")
    #plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig3.pdf", dpi=300, transparent=True)
    plt.show()