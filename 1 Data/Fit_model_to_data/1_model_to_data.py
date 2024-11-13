import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    parameters = [[10.2, 8.82e-03, 4.14e+02, 1.78e-01],
                 [14.1, 9.90e-03, 6.49e+02, 5.78e-02],
                 [11.6, 1.16e-02, 1.94e+03, 4.53e-02],
                 [5.34, 2.14e-02, 5.12e+02, 1.44e-01],
                 [4.20, 1.86e-02, 6.50e+02, 3.68e-02],
                 [12.9, 1.19e-02, 8.02e+02, 1.60e-01],
                 [2.71, 2.53e-02, 1.75e+02, 8.52e-02],
                 [5.48, 1.47e-02, 7.59e+02, 6.30e-02],
                 [5.48, 1.47e-02, 5.01e+02, 9.45e-02],
                 [5.48, 1.47e-02, 9.81e+02, 1.90e-02],
                 [52.1, 6.80e-03, 6.00e+02, 1.15e-01],
                 [6.46, 1.51e-02, 5.79e+02, 5.23e-02],
                 [11.4, 1.13e-02, 9.36e+02, 5.84e-02],
                 [2.43, 3.13e-02, 2.48e+02, 6.63e-02],
                 [12.2, 8.48e-03, 1.54e+02, 8.71e-02],
                 [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
                 [3.93, 2.05e-02, 1.60e+03, 4.10e-02],
                 [10.3, 8.65e-03, 4.15e+02, 5.27e-02],
                 [18.6, 6.02e-03, 9.31e+02, 4.84e-02],
                 [18.6, 6.02e-03, 9.40e+02, 4.83e-02],
                 [7.56, 1.17e-02, 3.65e+02, 8.35e-02],
                 [15.5, 7.47e-03, 1.92e+03, 3.41e-02],
                 [7.71, 1.09e-02, 5.37e+02, 2.38e-02],
                 [6.94, 1.74e-02, 1.19e+03, 5.05e-02]]

    data_idx = [1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                11, 12, 14,
                16, 19, 20,
                23, 27, 28,
                29, 30, 31,
                32, 33, 35]

    home = os.path.expanduser("~")
    for i in range(24):
        parameter = parameters[i]
        cell_idx = data_idx[i]

        st.set_default_plot_style()
        w  = 3.25*2.00
        h = 1.5*2.00
        fig = plt.figure(tight_layout=True, figsize=(w, h))
        gs = gridspec.GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharey=ax1)

        axis = [ax1, ax2]
        st.remove_top_right_axis(axis)
        ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
        ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')


        exp_isi = np.loadtxt(home + f"/Data/calcium/experimental/Spikes/HEK/HEK2/spike_times_{cell_idx:d}.dat")
        exp_idx = [i for i in range(len(exp_isi))]
        popt, pcov = curve_fit(fc.exponential_Ti, exp_idx, exp_isi, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        perr = np.sqrt(np.diag(pcov))
        mean_exp_isi = np.mean(exp_isi[int(2*popt[2]):])
        std_exp_isi = np.std(exp_isi[int(2*popt[2]):])
        cv_exp_isi = std_exp_isi / mean_exp_isi
        Tis = [popt[0] * np.exp(-idx / popt[2]) + popt[1] * (1 - np.exp(-idx / popt[2])) for idx in exp_idx]

        ax1.plot(exp_idx, Tis, c="k", lw=1,
                 label=f"$T_0$ = {popt[0]:.0f} $\pm$ {perr[0]:.0f}" + "\n" + f"$T_\infty$ = {popt[1]:.0f} $\pm$ {perr[1]:.0f}" + "\n" + rf"$n_{{\rm tr}}$ = {popt[2]:.1f} $\pm$ {perr[2]:.1f}")
        ax1.scatter(exp_idx, exp_isi, s=20, fc="w", ec=st.colors[5])
        ax1.axhline(popt[0], ls=":", c="C7")
        ax1.axhline(popt[1], ls=":", c="C7")
        ax1.axvspan(0, popt[2], facecolor="C7", alpha=0.5, zorder=0)

        ax1.set_ylabel("$T_{i}$")
        ax1.set_xlabel("$i$")
        ax1.set_ylim([0, 1.1 * max(exp_isi)])
        ax1.set_xlim([0, len(exp_isi)])
        legend = ax1.legend(loc=4, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 10})
        legend.get_frame().set_linewidth(0.5)

        tau = parameter[0]
        p = parameter[1]
        tau_er = parameter[2]
        eps_er = parameter[3]
        ax2.set_ylabel("$T_{i}$")
        ax2.set_xlabel("$i$")
        ax2.set_xlim([0, len(exp_isi)])
        data_sim_isi = np.loadtxt(home + f"/Data/calcium/markov/adaptive/fit/spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat")
        mean = np.mean(data_sim_isi)
        var = np.var(data_sim_isi)
        p1 = fc.k_corr(data_sim_isi, data_sim_isi, 1) / var
        data_sim_isi_tra = np.loadtxt(home + f"/Data/calcium/markov/adaptive/fit/transient_spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat")
        rows, cols = data_sim_isi_tra.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        mean_Ti = [np.mean(data_sim_isi_tra[:, idx]) for idx in idxs] # calculate the mean column wise

        for row in np.arange(1):
            ax2.scatter(idxs, data_sim_isi_tra[row, :], fc="w", ec=st.colors[0], alpha=0.75, s=20, zorder=1)
        ax2.scatter(idxs, mean_Ti, fc="w", ec="k", s=20, zorder=2, label=rf"CV = {np.sqrt(var)/mean:.2f}" + "\n" + rf"$\rho_1 = {p1:.2f}$")
        t0 = mean_Ti[0]
        t8 = mean_Ti[-1]
        func = lambda x, ntr: fc.exponential_Ti(x, t0, t8, ntr)  # fix t0, t8
        popt, pcov = curve_fit(func, idxs, mean_Ti, p0=[2.])
        ntr = popt[0]
        perr = np.sqrt(np.diag(pcov))
        ax2.axvspan(0, ntr, alpha=0.3, color="C7")

        # print(popt)
        ISI_fit = t0 * np.exp(-idxs / ntr) + t8 * (1. - np.exp(-idxs / ntr))
        ax2.plot(idxs, ISI_fit, c="k", zorder=3,
                 label=f"$T_0$ = {t0:.0f}" + "\n" + f"$T_\infty$ = {t8:.0f}" + "\n" + rf"$n_{{\rm tr}}$ = {ntr:.1f}")
        ax2.axhline(t0, ls=":", lw=1, c="k")
        ax2.axhline(t8, ls=":", lw=1, c="k")
        legend = ax2.legend(loc=4, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 10})
        legend.get_frame().set_linewidth(0.5)
        home = os.path.expanduser("~")
        plt.savefig(home + f"/Desktop/HEK Cells Fit/plot_fit/cell_{cell_idx:d}.png", dpi=300)
        plt.show()