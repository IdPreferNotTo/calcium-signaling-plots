import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit, least_squares

import styles as st
import functions as fc
import default_parameters as df

def exponential_cer(t, cer8, tau):
    return cer8 + (1 - cer8) * np.exp(-t / tau)

if __name__ == "__main__":
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.66 * 1.25
    fig, axs = plt.subplots(2, 2, layout="constrained", figsize=(w, h))
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[1, 0]
    ax4 = axs[1, 1]

    ax1.text(0.05, 0.95, r"A$_1$", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B$_1$", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"A$_2$", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"B$_2$", fontsize=11, transform=ax4.transAxes, va='top')

    axis_bot = [ax3, ax4]
    axis_top = [ax1, ax2]
    axis = axis_bot + axis_top
    st.remove_top_right_axis(axis)

    ax1.set_ylabel(r"$\Delta T$ / s")
    ax1.set_xlabel(r"$\tau_{er}$ / s")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax2.set_ylabel(r"$\Delta T$ / s")
    ax2.set_xlabel(r"$\varepsilon$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax3.set_ylabel(r"$n_{\rm tr}$")
    ax3.set_xlabel(r"$\tau_{er}$ / s")
    ax3.set_xscale("log")
    ax3.set_ylim([0, 5])

    ax4.set_ylabel(r"$n_{\rm tr}$")
    ax4.set_xlabel(r"$\varepsilon$")
    ax4.set_xscale("log")
    ax4.set_ylim([0, 7.5])

    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.03
    tau_er_fix = 300
    tau_ers = np.logspace(1, 3, 11)
    eps_ers = np.logspace(-2, 0, 11)

    home = os.path.expanduser("~")
    colors = [st.colors[1], st.colors[3]]
    for i, (tau, p), in enumerate(zip(taus, ps)):
        dTs = []
        n_trs = []
        dTs_langevin = []
        n_trs_langevin = []
        dTs_theory = []
        n_trs_theory = []
        n_trs_inspiration = []
        for tau_er in tau_ers:
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er_fix)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, bounds=([0, 0 ,0], [np.inf, np.inf, np.inf]))
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs.append(T8 - T0)
            n_trs.append(n_tr)

            data_isi = df.load_spike_times_langevin_transient(tau, p, tau_er, eps_er_fix)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, bounds=([0, 0 ,0], [np.inf, np.inf, np.inf]))
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs_langevin.append(T8 - T0)
            n_trs_langevin.append(n_tr)

            T0_theory = fc.calculate_T_init(tau, p)
            T8_theory = fc.self_consistent_T_infty(tau, p, tau_er, eps_er_fix)
            dTs_theory.append(T8_theory - T0_theory)
            tau_1, tau_2 = fc.calculate_tau_2(tau, p, tau_er, eps_er_fix)
            n_trs_theory.append(tau_1/T0_theory)
            n_trs_inspiration.append(tau_er/T8)

        if i == 0:
            label="mean-driven"
        else:
            label="excitable"
        ax1.scatter(tau_ers, dTs, fc="w", ec=colors[i], s=20, label=label)
        if i == 0:
            ax1.plot(tau_ers, dTs_theory, c=colors[i], ls=(0, (3, 1)))
        else:
            ax1.plot(tau_ers, dTs_theory, c=colors[i], ls=(0, (3, 1)), label="theory")

        ax3.scatter(tau_ers, n_trs, fc="w", ec=colors[i], s=20)


        ax3.plot(tau_ers, n_trs_theory, c=colors[i], zorder=3, ls=(0, (3, 1)))
        #ax3.plot(tau_ers, n_trs_inspiration, c=st.colors[5], zorder=3, ls=":")

        dTs = []
        n_trs = []
        dTs_langevin = []
        n_trs_langevin = []
        dTs_theory = []
        n_trs_theory = []
        n_trs_inspiration = []
        for eps_er in eps_ers:
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er_fix, eps_er)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, bounds=([0, 0 ,0], [np.inf, np.inf, np.inf]))
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs.append(T8 - T0)
            n_trs.append(n_tr)

            data_isi = df.load_spike_times_langevin_transient(tau, p, tau_er_fix, eps_er)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, bounds=([0, 0 ,0], [np.inf, np.inf, np.inf]))
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs_langevin.append(T8 - T0)
            n_trs_langevin.append(n_tr)

            T0_theory = fc.calculate_T_init(tau, p)
            T8_theory = fc.self_consistent_T_infty(tau, p, tau_er_fix, eps_er)
            dTs_theory.append(T8_theory - T0_theory)
            tau_1, tau_2 = fc.calculate_tau_2(tau, p, tau_er_fix, eps_er)
            n_trs_theory.append(tau_1/T0_theory)
            n_trs_inspiration.append(tau_er_fix/T8)

        ax4.scatter(eps_ers, n_trs, fc="w", ec=colors[i], s=20)
        ax2.scatter(eps_ers, dTs, fc="w", ec=colors[i], s=20)


        ax4.plot(eps_ers, n_trs_theory, c=colors[i], ls=(0, (3, 1)))
        #ax3.plot(eps_ers, n_trs_theory2, c=st.colors[5], ls="--")
        ax2.plot(eps_ers, dTs_theory, c=colors[i], ls=(0, (3, 1)))
        #ax4.plot(eps_ers, n_trs_inspiration, c=st.colors[5], zorder=3, ls=":")

    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.36, 0.0), loc=3,
                     ncol=4, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(1.)
    leg.legendHandles[2].set_color("k")

    home = os.path.expanduser("~")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig8.pdf", dpi=300, transparent=True)
    plt.show()

