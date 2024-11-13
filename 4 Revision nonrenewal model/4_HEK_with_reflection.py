import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

from scipy import stats
from scipy.optimize import curve_fit

import functions as fc
import styles as st

if __name__ == "__main__":
    home = os.path.expanduser("~")
    parameters_it0 = np.asarray([[10.2, 8.82e-03, 4.14e+02, 1.78e-01],
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
                 [6.94, 1.74e-02, 1.19e+03, 5.05e-02]])
    parameters_it2 = np.asarray([
        [6.09, 1.24e-02, 4.20e+02, 1.09e-01],
         [19.5, 8.95e-03, 7.06e+02, 6.44e-02],
         [22.4, 9.54e-03, 2.10e+03, 6.01e-02],
         [3.22, 2.73e-02, 5.18e+02, 9.78e-02],
         [7.11, 1.32e-02, 7.18e+02, 5.56e-02],
         [9.65, 1.31e-02, 8.00e+02, 1.31e-01],
         [1.84, 3.53e-02, 1.83e+02, 6.42e-02],
         [5.95, 1.40e-02, 7.87e+02, 6.71e-02],
         [4.42, 1.72e-02, 4.83e+02, 8.02e-02],
         [9.92, 1.04e-02, 1.04e+03, 3.13e-02],
         [53.3, 6.89e-03, 6.04e+02, 1.17e-01],
         [9.39, 1.26e-02, 5.91e+02, 6.72e-02],
         [23.5, 9.06e-03, 1.03e+03, 8.33e-02],
         [2.41, 3.14e-02, 2.28e+02, 6.75e-02],
         [10.2, 9.23e-03, 1.56e+02, 7.55e-02],
         [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
         [7.24, 1.42e-02, 1.69e+03, 6.79e-02],
         [13.9, 7.35e-03, 3.82e+02, 6.55e-02],
         [28.5, 5.03e-03, 9.54e+02, 5.86e-02],
         [35.4, 4.72e-03, 9.65e+02, 6.54e-02],
         [8.56, 1.08e-02, 3.53e+02, 9.42e-02],
         [32.3, 5.80e-03, 2.09e+03, 4.88e-02],
         [12.4, 8.30e-03, 5.46e+02, 3.36e-02],
         [77.0, 1.11e-02, 1.55e+03, 1.09e-01]])

    tau_ers = [p[2] for p in parameters_it2]
    eps_ers = [p[3] for p in parameters_it2]
    p1s_exp = []
    p1s_exp_err = []
    ntrs_exp = []
    dTs_exp = []

    stat_intervals = []
    p1s_sim = []
    ntrs_sim = []
    dTs_sim = []
    for tau, p, tau_er, eps_er in parameters_it2:
        file_seq = home + f"/Data/calcium/markov/adaptive/sequence_reflection/spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
        data_isi_seq = np.loadtxt(file_seq)
        covar = fc.k_corr(data_isi_seq[100:], data_isi_seq[100:], 1)
        var = fc.k_corr(data_isi_seq[100:], data_isi_seq[100:], 0)
        p1 = covar/var
        p1s_sim.append(p1)

        file2 = home +  f"/Data/calcium/markov/adaptive/transient_reflection/transient_spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
        data_isi_trans = np.loadtxt(file2)
        rows, cols = data_isi_trans.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        # calculate the mean column wise
        means_isi = np.mean(data_isi_trans, axis=0)
        t0 = means_isi[0]
        t8 = means_isi[-1]
        func = lambda x, ntr: fc.exponential_Ti(x, t0, t8, ntr)  # fix t0, t8
        popt, pcov = curve_fit(func, idxs, means_isi, p0=[2.])
        ntr = popt[0]
        dTs_sim.append(t8 - t0)
        ntrs_sim.append(ntr)

    data_idx = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 19, 20, 23, 27, 28, 29, 30, 31, 32, 33, 35]
    for i in data_idx:
        exp_isi = np.loadtxt(home + f"/Data/calcium/experimental/Spikes/HEK/HEK2/spike_times_{i:d}.dat")
        exp_idx = [i for i in range(len(exp_isi))]
        popt, pcov = curve_fit(fc.exponential_Ti, exp_idx, exp_isi, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        t0 = popt[0]
        t8 = popt[1]
        ntr = popt[2]
        dTs_exp.append(t8 - t0)
        ntrs_exp.append(ntr)
        idx_stat = int(2 * ntr)
        if len(exp_isi) < idx_stat + 1:
            continue
        stationary_isis = exp_isi[idx_stat:]
        var_0 = fc.k_corr(stationary_isis, stationary_isis, 0)
        var_1 = fc.k_corr(stationary_isis, stationary_isis, 1)
        p1 = var_1 / var_0
        print(i, len(stationary_isis), f"{p1:.2f}")
        stat_intervals.append(len(stationary_isis))
        p1s_exp.append(p1)
        p1s_exp_err.append(np.sqrt(1 / (len(exp_isi[idx_stat:]) - 1) * (1. + np.power(p1, 2))))

    print(np.mean(stat_intervals), np.std(stat_intervals))
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 3.75 * 1.25
    # fig = plt.figure(tight_layout=True, figsize=(w, h))
    # gs = gridspec.GridSpec(2, 2)
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax2 = fig.add_subplot(gs[0, 1])
    # ax3 = fig.add_subplot(gs[1, 0])
    # ax4 = fig.add_subplot(gs[1, 1])
    fig, axs = plt.subplots(3, 2, layout="constrained", figsize=(w, h))
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[1, 0]
    ax4 = axs[1, 1]
    ax5 = axs[2, 0]
    ax6 = axs[2, 1]
    axis = [ax1, ax2, ax3, ax4, ax5, ax6]
    st.remove_top_right_axis(axis)
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"C", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"D", fontsize=11, transform=ax4.transAxes, va='top')
    ax5.text(0.05, 0.95, r"E", fontsize=11, transform=ax5.transAxes, va='top')
    ax6.text(0.05, 0.95, r"F", fontsize=11, transform=ax6.transAxes, va='top')
    # ax1.scatter(tau_ers, ntrs, s=20, fc="w", ec=st.colors[0])
    # res_p = stats.pearsonr(tau_ers, ntrs)
    # print(res_p)
    ax1.set_xlabel(r"$\Delta T$ / s")
    ax1.set_ylabel(r"$\rho_1$")
    ax1.set_ylim([-0.5, 0.15])
    ax1.scatter(dTs_sim, p1s_sim, s=20, fc="w", ec=st.colors[0])
    res_l = stats.linregress(dTs_sim, p1s_sim)
    res_p = stats.pearsonr(dTs_sim, p1s_sim)
    pcor = res_p.statistic
    xs = np.linspace(0., 400, 100)
    ax1.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"Pearson = ${pcor:.1f}$")
    ax5.plot(xs, res_l.intercept + res_l.slope * xs, color=st.colors[0])
    ax1.axhline(0, ls=":", c="k")
    leg1 = ax1.legend(fancybox=False, edgecolor="k", framealpha=1, fontsize=7, loc=1)
    leg1.get_frame().set_linewidth(0.75)

    ax2.set_xlabel(r"$n_{\rm tr}$")
    ax2.set_ylabel(r"$\rho_1$")
    ax2.set_ylim([-0.5, 0.15])
    ax2.scatter(ntrs_sim, p1s_sim, s=20, fc="w", ec=st.colors[0])
    ntrs_sim = np.asarray(ntrs_sim)
    res_l = stats.linregress(ntrs_sim, p1s_sim)
    res_p = stats.pearsonr(ntrs_sim, p1s_sim)
    pcor = res_p.statistic
    xs = np.linspace(0., 10, 100)
    ax2.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"Pearson $= {pcor:.1f}$")
    ax6.plot(xs, res_l.intercept + res_l.slope * xs, color=st.colors[0])
    ax2.axhline(0, ls=":", c="k")
    leg2 = ax2.legend(fancybox=False, edgecolor="k", framealpha=1, fontsize=7, loc=1)
    leg2.get_frame().set_linewidth(0.75)

    ax3.set_xlabel(r"$\Delta T$ / s")
    ax3.set_ylabel(r"$\rho_1$")
    ax3.set_ylim([-0.75, 0.75])
    ax3.scatter(dTs_exp, p1s_exp, s=20, fc="w", ec=st.colors[7], zorder=2)
    for dT, p1, err in zip(dTs_exp, p1s_exp, p1s_exp_err):
        ax3.plot([dT, dT], [p1 - err, p1 + err], lw=1, color=st.colors[7], zorder=1)
    res_l = stats.linregress(dTs_exp, p1s_exp)
    res_p = stats.pearsonr(dTs_exp, p1s_exp)
    pcor = res_p.statistic
    xs = np.linspace(0., 400, 100)
    ax3.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"Pearson = ${pcor:.1f}$")
    ax5.plot(xs, res_l.intercept + res_l.slope * xs, color=st.colors[7], label=rf"Pearson = ${pcor:.1f}$")
    ax3.axhline(0, ls=":", c="k", zorder=1)
    leg3 = ax3.legend(fancybox=False, edgecolor="k", framealpha=1, fontsize=7, loc=1)
    leg3.get_frame().set_linewidth(0.75)

    ax4.set_xlabel(r"$n_{\rm tr}$")
    ax4.set_ylabel(r"$\rho_1$")
    ax4.set_ylim([-0.75, 0.75])
    ax4.scatter(ntrs_exp, p1s_exp, s=20, fc="w", ec=st.colors[7], zorder=2)
    for ntr, p1, err in zip(ntrs_exp, p1s_exp, p1s_exp_err):
        ax4.plot([ntr, ntr], [p1 - err, p1 + err], lw=1, color=st.colors[7], zorder=1)
    res_l = stats.linregress(ntrs_exp, p1s_exp)
    res_p = stats.pearsonr(ntrs_exp, p1s_exp)
    pcor = res_p.statistic
    xs = np.linspace(0., 10, 100)
    ax4.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"Pearson = ${pcor:.1f}$")
    ax6.plot(xs, res_l.intercept + res_l.slope * xs, color=st.colors[7])
    ax4.axhline(0, ls=":", c="k", zorder=1)
    leg4 = ax4.legend(fancybox=False, edgecolor="k", framealpha=1, fontsize=7, loc=1)
    leg4.get_frame().set_linewidth(0.75)


    #plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/SUB2/figures/fig12.pdf", dpi=300, transparent=True)
    plt.show()
