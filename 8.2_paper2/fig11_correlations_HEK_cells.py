import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.optimize import curve_fit

import functions as fc
import styles as st


def transient_func(i, T0, T8, iota):
    return T0 * np.exp(-i / iota) + T8 * (1. - np.exp(-i / iota))


if __name__ == "__main__":
    home = os.path.expanduser("~")

    ntrs = []
    dTs = []
    p1s = []
    err_p1s = []
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
        print(i, p1)
        ntrs.append(ntr)
        dTs.append(t8 - t0)
        p1s.append(p1)
        err_p1s.append(np.sqrt(1 / (len(Tis[idx_stat:]) - 1) * (1. + np.power(p1, 2))))

    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 1.5 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    axis = [ax1, ax2]
    st.remove_top_right_axis(axis)
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')

    # ax1.scatter(tau_ers, ntrs, s=20, fc="w", ec=st.colors[0])
    # res_p = stats.pearsonr(tau_ers, ntrs)
    # print(res_p)
    ax1.set_xlabel(r"$\Delta T$")
    ax1.set_ylabel(r"$\rho_1$")
    ax1.set_ylim([-1.0, 1.0])
    ax1.scatter(dTs, p1s, ec=st.colors[5], fc="w", s=15, zorder=2)
    for dT, p1, err in zip(dTs, p1s, err_p1s):
        ax1.plot([dT, dT], [p1 - err, p1 + err], lw=1, color=st.colors[5], zorder=1)
    res_l = stats.linregress(dTs, p1s)
    res_p = stats.pearsonr(dTs, p1s)
    print(res_p)
    pcor = res_p.statistic
    xs = np.linspace(0., 400, 100)
    ax1.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax1.axhline(0, ls=":", c="k")
    leg1 = ax1.legend(fancybox=False, edgecolor="k", fontsize=7, framealpha=1.)
    leg1.get_frame().set_linewidth(0.75)

    ax2.set_xlabel(r"$n_{\rm tr}$")
    ax2.set_ylabel(r"$\rho_1$")
    ax2.set_ylim([-1.0, 1.0])
    ax2.scatter(ntrs, p1s, ec=st.colors[5], fc="w", s=15, zorder=2)
    for ntr, p1, err in zip(ntrs, p1s, err_p1s):
        ax2.plot([ntr, ntr], [p1 - err, p1 + err], lw=1, color=st.colors[5], zorder=1)

    ntrs = np.asarray(ntrs)
    res_l = stats.linregress(ntrs, p1s)
    res_p = stats.pearsonr(ntrs, p1s)
    print(res_p)
    pcor = res_p.statistic
    xs = np.linspace(0., 10, 100)
    ax2.plot(xs, res_l.intercept + res_l.slope * xs, 'r', label=rf"$\rho = {pcor:.1f}$")
    ax2.axhline(0, ls=":", c="k")
    leg2 = ax2.legend(fancybox=False, edgecolor="k", fontsize=7, framealpha=1.)
    leg2.get_frame().set_linewidth(0.75)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig11.pdf", transparent=True)
    plt.show()
