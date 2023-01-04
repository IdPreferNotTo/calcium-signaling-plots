import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 3))
    gs = gridspec.GridSpec(2, 3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])
    axis_top = [ax1, ax2, ax3]
    axis_bot = [ax4, ax5, ax6]
    axis = axis_top + axis_bot
    st.remove_top_right_axis(axis)

    ax1.set_ylabel(r"$\langle T \rangle$")
    ax1.set_xlabel(r"$\tau_{er}$")
    #ax1.set_xscale("log")

    ax2.set_ylabel(r"$CV_T$")
    ax2.set_xlabel(r"$\tau_{er}$")
    #ax2.set_xscale("log")

    ax3.set_ylabel(r"$\rho_1$")
    ax3.set_xlabel(r"$\tau_{er}$")

    ax4.set_ylabel(r"$\langle T \rangle$")
    ax4.set_xlabel(r"$\varepsilon_{er}$")
    #ax3.set_xscale("log")

    ax5.set_ylabel(r"$CV_T$")
    ax5.set_xlabel(r"$\varepsilon_{er}$")
    #ax4.set_xscale("log")

    ax6.set_ylabel(r"$\rho_1$")
    ax6.set_xlabel(r"$\varepsilon_{er}$")

    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.1
    tau_er_fix = 100
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)

    for i, (tau, p) in enumerate(zip(taus, ps)):
        means = []
        cvs = []
        p1s = []
        for tau_er in tau_ers:
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er_fix)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            cvs.append(cv)
            p1s.append(p1)

        ax1.scatter(tau_ers, means, fc="w", ec=st.colors[i], s=20, zorder=3)
        ax2.scatter(tau_ers, cvs, fc="w", ec=st.colors[i], s=20, zorder=3)
        ax3.scatter(tau_ers, p1s, fc="w", ec=st.colors[i], s=20, zorder=3)

        means = []
        cvs = []
        p1s = []
        for eps_er in eps_ers:
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er_fix, ampa=eps_er)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            cvs.append(cv)
            p1s.append(p1)
        ax4.scatter(eps_ers, means, fc="w", ec=st.colors[i], s=20, zorder=3)
        ax5.scatter(eps_ers, cvs, fc="w", ec=st.colors[i], s=20, zorder=3)
        ax6.scatter(eps_ers, p1s, fc="w", ec=st.colors[i], s=20, zorder=3)

    home = os.path.expanduser("~")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig5.pdf", transparent=True)
    plt.show()