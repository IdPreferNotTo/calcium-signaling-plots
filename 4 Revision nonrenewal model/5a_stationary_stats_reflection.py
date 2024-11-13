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
    w = 2.0 * 1.25
    h = 4.0 * 1.25
    # fig = plt.figure(tight_layout=True, figsize=(w, h))
    # gs = gridspec.GridSpec(3, 2)
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax2 = fig.add_subplot(gs[1, 0])
    # ax3 = fig.add_subplot(gs[2, 0])
    # ax4 = fig.add_subplot(gs[0, 1])
    # ax5 = fig.add_subplot(gs[1, 1])
    # ax6 = fig.add_subplot(gs[2, 1])

    fig, axs = plt.subplots(nrows=3, ncols=1, layout="constrained", figsize=(w, h))
    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]

    axis = [ax1, ax2, ax3]
    st.remove_top_right_axis(axis)
    ax1.text(0.15, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.15, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.15, 0.95, r"C", fontsize=11, transform=ax3.transAxes, va='top')



    ax1.set_ylabel(r"$\langle T \rangle$ / s")
    ax1.set_xlabel(r"$\varepsilon$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylim([20, 2000])

    ax2.set_ylabel(r"$C_V$")
    ax2.set_xlabel(r"$\varepsilon$")
    ax2.set_xscale("log")
    ax2.set_ylim([0, 0.8])

    ax3.set_ylabel(r"$\rho_1$")
    ax3.set_xlabel(r"$\varepsilon$")
    ax3.set_xscale("log")
    ax3.set_ylim([-0.5, 0.2])
    ax3.axhline(0, ls=":", c="k")

    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    tau_ers = np.logspace(1, 3, 11)
    eps_ers = np.logspace(-2, 0, 11)

    #colors = ["#007A45", "#20639B"]
    colors = [st.colors[1], st.colors[3]]
    for i, (tau, p) in enumerate(zip(taus, ps)):
        if i == 0:
            label = "mean-driven"
        else:
            label = "excitable"

        means_d0 = []
        cvs_d0 = []
        p1s_d0 = []
        for ii, eps_er in enumerate(eps_ers):
            data_isi_d0 = df.load_spike_times_markov(tau, p, cer=True, taua=300, ampa=eps_er)
            mean = np.mean(data_isi_d0)
            std = np.std(data_isi_d0)
            var = np.var(data_isi_d0)
            cv = std / mean
            p1 = fc.k_corr(data_isi_d0, data_isi_d0, 1) / var
            means_d0.append(mean)
            cvs_d0.append(cv)
            p1s_d0.append(p1)

        if i == 0:
            axis[0].plot(eps_ers, means_d0, c=colors[i], zorder=1, label=f"$\delta = 0$")
        else:
            axis[0].plot(eps_ers, means_d0, c=colors[i], zorder=1)
        axis[1].plot(eps_ers, cvs_d0, c=colors[i], zorder=1)
        axis[2].plot(eps_ers, p1s_d0, c=colors[i], zorder=1)

        for n, d in enumerate([0.01]):
            means = []
            cvs = []
            p1s = []
            for ii, eps_er in enumerate(eps_ers):
                data_isi = np.loadtxt(home + f"/Data/calcium/markov/adaptive/sequence_reflection/spike_times_markov_ip1.00_taua3.00e+02_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5_d{d}.dat")
                mean = np.mean(data_isi)
                std = np.std(data_isi)
                var = np.var(data_isi)
                cv = std/mean
                p1 = fc.k_corr(data_isi, data_isi, 1)/var
                means.append(mean)
                cvs.append(cv)
                p1s.append(p1)

            if i == 0:
                axis[0].plot(eps_ers, means, c=colors[i], alpha= 0.3, zorder=2, label=f"$\delta = {d}$")
            else:
                axis[0].plot(eps_ers, means, c=colors[i], alpha=0.3, zorder=2)
            axis[1].plot(eps_ers, cvs, c=colors[i], alpha=0.3, zorder=2)
            axis[2].plot(eps_ers, p1s, c=colors[i], alpha=0.3, zorder=2)



    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 1.0, .0), loc=3,
                      ncol=4, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(.75)
    #leg.legendHandles[2].set_color("k")
    home = os.path.expanduser("~")
    plt.savefig(home + f"/Desktop/hek_fit/stationary_statistics_with_reflection.pdf", dpi=300)
    plt.show()