import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from functions import *


def gauss_dist(xs, mean, std):
    gauss_dist = []
    for x in xs:
        gauss = 1 / np.sqrt(2 * np.pi * (std ** 2)) * np.exp(-((x - mean) ** 2) / (2 * std ** 2))
        gauss_dist.append(gauss)
    return gauss_dist

if __name__ == "__main__":
    # Set up Data
    ca_fix = 0.50
    ca_res = 0.33
    n_cl = 10
    n_ch = 4
    jca = 1
    tau = 1
    home = os.path.expanduser("~")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/"
    file = "ca_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca_fix)
    data = np.loadtxt(home + folder + file)
    ts, cas, jpuffs, adaps = np.transpose(data)
    ts_clear = []
    cas_clear = []
    jpuffs_clear = []
    mean_jpuffs_clear = []
    for t, ca, jpuff in zip(ts, cas, jpuffs):
        if ca != 1:
            ts_clear.append(t)
            cas_clear.append(ca)
            jpuffs_clear.append(jpuff)
            Popen = np.power(1 + 2 * np.power(n_ch * (n_ch + 1), -1, dtype=float) * np.power(ca_res / ca, 3) * (1 + ca ** 3) / (1 + ca_res ** 3) * 10 / 0.02, -1)
            mean_jpuff = jca * (n_cl * (n_ch + 2) / 3) * Popen
            mean_jpuffs_clear.append(mean_jpuff)
    jpuffs_dt01 = jpuffs_clear
    jpuffs_dt1 = coarse_grain_list(jpuffs_clear, 10)

    # Set up plot style
    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(2, 2)
    ax_ts1 = fig.add_subplot(gs[0, 0])
    ax_hist1 = fig.add_subplot(gs[0, 1])
    ax_ts2 = fig.add_subplot(gs[1, 0])
    ax_hist2 = fig.add_subplot(gs[1, 1])
    remove_top_right_axis([ax_ts1, ax_hist1, ax_ts2, ax_hist2])

    plt_start = 0
    plt_stop = plt_start + 500

    ax_ts1.set_xlabel("$t$")
    ax_ts1.set_xlim([0, 50])
    ax_ts1.set_ylabel(r"$y_{0.1}$")
    ax_ts1.set_ylim([0, 6])
    mean = np.mean(jpuffs_clear)
    ax_ts1.plot(ts_clear[plt_start:plt_stop], jpuffs_clear[plt_start:plt_stop], c="C3")
    ax_ts1.plot(ts_clear[plt_start:plt_stop], mean_jpuffs_clear[plt_start:plt_stop], c="k",
                label="$\mu = {:.2f}$".format(mean))
    legend_ts1 = ax_ts1.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_ts1.get_frame().set_linewidth(0.5)

    ts_clear_dt1 = [10 * t for t in ts_clear]
    ts_clear_dt1_new = []
    jpuffs_dt1_new = []
    for jpuff1, t1, jpuff2, t2 in zip(jpuffs_dt1[:-1], ts_clear_dt1[:-1], jpuffs_dt1[1:],
                                      ts_clear_dt1[1:]):  # set  = [time, state, idx]
        ts_clear_dt1_new.append(t1)
        jpuffs_dt1_new.append(jpuff1)
        ts_clear_dt1_new.append(t2)
        jpuffs_dt1_new.append(jpuff1)
    mean = np.mean(jpuffs_clear)

    ax_ts2.set_xlabel("$t$")
    ax_ts2.set_xlim([0, 50])
    ax_ts2.set_ylabel(r"$y_{1.0}$")
    ax_ts2.set_ylim([0, 6])


    ax_ts2.plot(ts_clear_dt1_new[plt_start:plt_stop], jpuffs_dt1_new[plt_start:plt_stop], c="C0")
    ax_ts2.plot(ts_clear_dt1[plt_start:plt_stop], mean_jpuffs_clear[plt_start:plt_stop], c="k",
                label="$\mu = {:.2f}$".format(mean))

    legend_ts2 = ax_ts2.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_ts2.get_frame().set_linewidth(0.5)

    std01 = np.std(jpuffs_dt01)
    std1 = np.std(jpuffs_dt1)
    var01 = np.var(jpuffs_dt01)
    var1 = np.var(jpuffs_dt1)
    gauss001 = np.linspace(mean - 3 * std01, mean + 3 * std01, 100)
    gauss01 = np.linspace(mean - 3 * std1, mean + 3 * std1, 100)

    ax_hist1.set_xlabel("$y_{0.1}$")
    ax_hist1.set_xlim([-1, 6])
    ax_hist1.set_ylabel("$P(y_{0.1})$")
    ax_hist1.set_ylim([0, 1.25])
    ax_hist1.plot(gauss001, gauss_dist(gauss001, mean, std01), c="C3")
    ax_hist1.hist(jpuffs_dt01, bins=50, alpha=0.7, color="C3", density=True,
                  label=r"$\Delta t = {:.1f}$".format(0.1) + "\n" + "$\sigma_{y_{0.1}}$" + " = {:.2f}".format(std01))
    legend_hs1 = ax_hist1.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_hs1.get_frame().set_linewidth(0.5)


    ax_hist2.set_xlabel("$y_{1.0}$")
    ax_hist2.set_xlim([-1, 6])
    ax_hist2.set_ylabel("$P(y_{1.0})$")
    ax_hist2.set_ylim([0, 1.25])
    ax_hist2.plot(gauss01, gauss_dist(gauss01, mean, std1), c="C0")
    ax_hist2.hist(jpuffs_dt1, bins=50, alpha=0.7, color="C0", density=True,
                  label=r"$\Delta t = {:.1f}$".format(1.0) + "\n" + "$\sigma_{y_{1.0}}$" + " = {:.2f}".format(std1))
    legend_hs2 = ax_hist2.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_hs2.get_frame().set_linewidth(0.5)

    plt.savefig(home + "/Data/Calcium/Plots/markov_jpuff_coarse_grained_static.pdf".format(ca_fix, tau, jca),
                transparent=True)
    plt.show()
