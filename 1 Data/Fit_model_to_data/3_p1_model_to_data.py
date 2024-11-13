import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from scipy import stats

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    home = os.path.expanduser("~")
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

    p1s = []
    ntrs = []
    dTs = []
    for tau, p, tau_er, eps_er in parameters:
        file1 = home + f"/Data/calcium/markov/adaptive/fit/spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
        data_isi = np.loadtxt(file1)
        covar = fc.k_corr(data_isi[100:], data_isi[100:], 1)
        var = fc.k_corr(data_isi[100:], data_isi[100:], 0)
        p1 = covar/var
        p1s.append(p1)

        file2 = home + f"/Data/calcium/markov/adaptive/fit/transient_spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
        data_isi_trans = np.loadtxt(file2)
        rows, cols = data_isi_trans.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        means_isi = [np.mean(data_isi_trans[:, idx]) for idx in idxs] # calculate the mean column wise

        popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_isi, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        dTs.append(popt[1] - popt[0])
        ntrs.append(popt[2])

    st.set_default_plot_style()
    w  = 3.25*1.25
    h = 1.5*1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    axis = [ax1, ax2]
    st.remove_top_right_axis(axis)
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')

    ax1.set_xlabel(r"$\Delta T$")
    ax1.set_ylabel(r"$\rho_1$")
    # ax1.set_xlim([0., 0.5])
    ax1.set_ylim([-0.5, 0.1])
    ax1.scatter(dTs, p1s, s=20, fc="w", ec=st.colors[0])
    res_l = stats.linregress(dTs, p1s)
    res_p = stats.pearsonr(dTs, p1s)
    pcor = res_p.statistic
    xs = np.linspace(0., 400, 100)
    ax1.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax1.axhline(0, ls=":", c="k")
    ax1.legend(fancybox=False, fontsize=8)

    ax2.set_xlabel(r"$n_{\rm tr}$")
    ax2.set_ylabel(r"$\rho_1$")
    ax2.set_ylim([-0.5, 0.1])
    ax2.scatter(ntrs, p1s, s=20, fc="w", ec=st.colors[0])
    ntrs = np.asarray(ntrs)
    res_l = stats.linregress(ntrs, p1s)
    res_p = stats.pearsonr(ntrs, p1s)
    pcor = res_p.statistic
    xs = np.linspace(0., 10, 100)
    ax2.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax2.axhline(0, ls=":", c="k")
    ax2.legend(fancybox=False, fontsize=8)

    plt.show()
