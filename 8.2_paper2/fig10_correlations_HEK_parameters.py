import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.optimize import curve_fit

import functions as fc
import styles as st

if __name__ == "__main__":
    home = os.path.expanduser("~")
    parameter2 = [[1.01e+01, 8.46e-03, 4.79e+02, 2.83e-01],
                  [1.25e+01, 9.07e-03, 6.69e+02, 1.05e-01],
                  [7.50e+00, 1.19e-02, 1.41e+03, 5.84e-02],
                  [5.42e+00, 1.41e-02, 5.06e+02, 9.97e-02],
                  [4.82e+00, 1.54e-02, 6.88e+02, 5.27e-02],
                  [1.02e+01, 9.32e-03, 8.46e+02, 1.68e-01],
                  [3.04e+00, 2.24e-02, 1.91e+02, 1.31e-01],
                  [5.74e+00, 1.35e-02, 7.59e+02, 1.06e-01],
                  [5.74e+00, 1.35e-02, 5.03e+02, 1.60e-01],
                  [6.68e+00, 1.17e-02, 8.27e+02, 3.53e-02],
                  [7.32e+01, 5.05e-03, 9.10e+02, 2.43e-01],
                  [6.90e+00, 1.19e-02, 5.70e+02, 6.29e-02],
                  [8.72e+00, 1.02e-02, 8.46e+02, 6.76e-02],
                  [3.45e+00, 2.01e-02, 2.38e+02, 6.10e-02],
                  [1.20e+01, 8.53e-03, 1.76e+02, 2.01e-01],
                  [3.60e+00, 1.94e-02, 3.04e+02, 7.38e-02],
                  [5.66e+00, 1.32e-02, 1.61e+03, 5.26e-02],
                  [1.01e+01, 8.46e-03, 4.04e+02, 1.00e-01],
                  [1.97e+01, 6.36e-03, 9.73e+02, 1.48e-01],
                  [1.97e+01, 6.36e-03, 9.46e+02, 1.50e-01],
                  [7.94e+00, 1.03e-02, 3.53e+02, 1.41e-01],
                  [1.57e+01, 7.55e-03, 1.63e+03, 9.84e-02],
                  [7.94e+00, 1.03e-02, 5.16e+02, 4.67e-02],
                  [7.19e+00, 1.11e-02, 1.03e+03, 3.96e-02]] #tau, p, tau_er, eps
    statistics2 = np.asarray([[45, 0.12, 1.2, 343, -0.38],
                       [30, 0.11, 5.0, 122, 0.17],
                       [30, 0.14, 8.1, 228, -0.10],
                       [35, 0.18, 2.3, 219, -0.38],
                       [35, 0.20, 3.8, 188, -0.22],
                       [35, 0.12, 2.6, 285, -0.39],
                       [35, 0.31, 1.0, 187, 0.24],
                       [35, 0.17, 2.4, 314, -0.38],
                       [35, 0.17, 1.6, 297, -0.28],
                       [40, 0.17, 7.3, 148, 0.19],
                       [35, 0.07, 3.6, 139, 0.49],
                       [35, 0.15, 4.4, 137, -0.29],
                       [35, 0.13, 5.6, 165, 0.11],
                       [35, 0.27, 2.0, 123, -0.19],
                       [35, 0.11, 1.7, 91, 0.30],
                       [35, 0.26, 1.9, 165, -0.42],
                       [40, 0.18, 4.3, 414, -0.61],
                       [45, 0.12, 2.7, 143, 0.27],
                       [40, 0.09, 4.1, 193, 0.27],
                       [40, 0.09, 4.0, 193, 0.19],
                       [40, 0.14, 1.9, 174, 0.18],
                       [35, 0.10, 7.4, 232, 0.42],
                       [40, 0.14, 5.1, 110, -0.34],
                       [40, 0.15, 7.5, 181, -0.23]])     #t0, cv, ntr, t8, p1

    tau_ers = [p[2] for p in parameter2]
    eps_ers = [p[3] for p in parameter2]
    p1s_exp = []
    p1s_exp_err = []
    p1s = []
    ntrs = []
    dTs = []
    for tau, p, tau_er, eps_er in parameter2:
        file1 = home + f"/Desktop/HEK Cells Fit/data_fix_t0/spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
        data_isi = np.loadtxt(file1)
        covar = fc.k_corr(data_isi[100:], data_isi[100:], 1)
        var = fc.k_corr(data_isi[100:], data_isi[100:], 0)
        p1 = covar/var
        p1s.append(p1)

        file2 = home + f"/Desktop/HEK Cells Fit/data_fix_t0/transient_spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
        data_isi_trans = np.loadtxt(file2)
        rows, cols = data_isi_trans.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        means_isi = [np.mean(data_isi_trans[:, idx]) for idx in idxs] # calculate the mean column wise

        t0 = means_isi[0]
        func = lambda x, t8, ntr: fc.exponential_Ti(x, t0, t8, ntr)
        popt, pcov = curve_fit(func, idxs, means_isi, p0=(means_isi[-1], 5.))
        dTs.append(popt[0] - t0)
        ntrs.append(popt[1])

    for i in range(1, 39):
        if i in [10, 13, 15, 17, 18, 21, 22, 24, 25, 26, 34, 36, 37, 38]:
            continue

        file = home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/spike_times_{i:d}.dat"
        Tis = np.loadtxt(file)  # interspike intervals
        t0 = Tis[0]
        t8 = Tis[-1]
        idxs = np.arange(len(Tis))
        func = lambda x, ntr, t8: fc.exponential_Ti(x, t0, t8, ntr)  # fix t0, t8
        popt, pcov = curve_fit(func, idxs, Tis, p0=[t8, 2.])
        ntr = popt[0]
        t8 = popt[1]

        idx_stat = int(2 * ntr)
        if len(Tis) < idx_stat + 1:
            continue

        stationary_isis = Tis[idx_stat:]
        var_0 = fc.k_corr(stationary_isis, stationary_isis, 0)
        var_1 = fc.k_corr(stationary_isis, stationary_isis, 1)
        p1 = var_1 / var_0
        p1s_exp.append(p1)
        p1s_exp_err.append(np.sqrt(1 / (len(Tis[idx_stat:]) - 1) * (1. + np.power(p1, 2))))

    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.66 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axis = [ax1, ax2, ax3, ax4]
    st.remove_top_right_axis(axis)
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"C", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"D", fontsize=11, transform=ax4.transAxes, va='top')
    # ax1.scatter(tau_ers, ntrs, s=20, fc="w", ec=st.colors[0])
    # res_p = stats.pearsonr(tau_ers, ntrs)
    # print(res_p)
    ax1.set_xlabel(r"$\Delta T$")
    ax1.set_ylabel(r"$\rho_1$")
    ax1.set_ylim([-0.5, 0.15])
    ax1.scatter(dTs, p1s, s=20, fc="w", ec=st.colors[0])
    res_l = stats.linregress(dTs, p1s)
    res_p = stats.pearsonr(dTs, p1s)
    pcor = res_p.statistic
    xs = np.linspace(0., 400, 100)
    ax1.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax1.axhline(0, ls=":", c="k")
    leg1 = ax1.legend(fancybox=False, edgecolor="k", fontsize=7)
    leg1.get_frame().set_linewidth(0.75)

    ax2.set_xlabel(r"$n_{\rm tr}$")
    ax2.set_ylabel(r"$\rho_1$")
    ax2.set_ylim([-0.5, 0.15])
    ax2.scatter(ntrs, p1s, s=20, fc="w", ec=st.colors[0])
    ntrs = np.asarray(ntrs)
    res_l = stats.linregress(ntrs, p1s)
    res_p = stats.pearsonr(ntrs, p1s)
    pcor = res_p.statistic
    xs = np.linspace(0., 10, 100)
    ax2.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax2.axhline(0, ls=":", c="k")
    leg2 = ax2.legend(fancybox=False, edgecolor="k", fontsize=7)
    leg2.get_frame().set_linewidth(0.75)

    ax3.set_xlabel(r"$\Delta T$")
    ax3.set_ylabel(r"$\rho_1$")
    ax3.scatter(dTs, p1s_exp, s=20, fc="w", ec=st.colors[5], zorder=2)
    for dT, p1, err in zip(dTs, p1s_exp, p1s_exp_err):
        ax3.plot([dT, dT], [p1 - err, p1 + err], lw=1, color=st.colors[5], zorder=1)
    res_l = stats.linregress(dTs, p1s_exp)
    res_p = stats.pearsonr(dTs, p1s_exp)
    pcor = res_p.statistic
    xs = np.linspace(0., 400, 100)
    ax3.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax3.axhline(0, ls=":", c="k")
    leg3 = ax3.legend(fancybox=False, edgecolor="k", fontsize=7)
    leg3.get_frame().set_linewidth(0.75)

    ax4.set_xlabel(r"$n_{\rm tr}$")
    ax4.set_ylabel(r"$\rho_1$")
    ax4.scatter(ntrs, p1s_exp, s=20, fc="w", ec=st.colors[5], zorder=2)
    for ntr, p1, err in zip(ntrs, p1s_exp, p1s_exp_err):
        ax4.plot([ntr, ntr], [p1 - err, p1 + err], lw=1, color=st.colors[5], zorder=1)
    ntrs = np.asarray(ntrs)
    res_l = stats.linregress(ntrs, p1s_exp)
    res_p = stats.pearsonr(ntrs, p1s_exp)
    pcor = res_p.statistic
    xs = np.linspace(0., 10, 100)
    ax4.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax4.axhline(0, ls=":", c="k")
    leg4 = ax4.legend(fancybox=False, edgecolor="k", fontsize=7)
    leg4.get_frame().set_linewidth(0.75)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig10.pdf", transparent=True)
    plt.show()
