import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    # Set up Parameters
    ca = 0.50
    ca_rest = 0.33
    N = 10
    n = 5
    m = 4
    jca = 1
    tau = 1

    # Load Data
    home = os.path.expanduser("~")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/ca_fix/"
    file = f"ca_markov_cafix{ca:.2f}_ip1.00_tau{tau:.2e}_j{jca:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    ts, cas, jpuffs, adaps = np.transpose(data)

    r_opn = 0.13 * np.power(ca / ca_rest, 3) * (1 + ca_rest ** 3) / (1 + ca ** 3)
    r_ref = 1.3 * np.power(ca / ca_rest, 3) * (1 + ca_rest ** 3) / (1 + ca ** 3)
    r_cls = 50
    p0s = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)
    xs = np.empty(n + m)
    for i in range(n + m):
        if i < m:
            xs[i] = 0
        if i >= m:
            xs[i] = n + m - i
    mean_x = np.dot(xs, p0s)

    # alternatively one could use an analytical expression for the mean puff current
    Popen = np.power(1 + 2 * np.power(n * (n + 1), -1, dtype=float) * np.power(ca_rest / ca, 3) * (1 + ca ** 3) / (1 + ca_rest ** 3) * 10 / 0.02, -1)
    mean_jpuff = jca * (N * (n + 2) / 3) * Popen

    jpuffs_dt01 = jpuffs
    jpuffs_dt1 = fc.coarse_grain_list(jpuffs, 10)

    # Set up plot style
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(2, 2)
    ax_ts1 = fig.add_subplot(gs[0, 0])
    ax_hist1 = fig.add_subplot(gs[0, 1])
    ax_ts2 = fig.add_subplot(gs[1, 0])
    ax_hist2 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax_ts1, ax_hist1, ax_ts2, ax_hist2])

    plt_start = 0
    plt_stop = plt_start + 500

    ax_ts1.set_xlabel("$t$")
    ax_ts1.set_xlim([0, 50])
    ax_ts1.set_ylabel(r"$y_{0.1}$")
    ax_ts1.set_ylim([0, 6])
    mean = np.mean(jpuffs)
    ax_ts1.plot(ts[plt_start:plt_stop], jpuffs[plt_start:plt_stop], c="C3")
    ax_ts1.axhline(N*mean_x, c="k", label="$\mu = {:.2f}$".format(mean))
    legend_ts1 = ax_ts1.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_ts1.get_frame().set_linewidth(0.5)

    ts_dt1 = [10 * t for t in ts]
    ts_dt1_new = []
    jpuffs_dt1_new = []
    for jpuff1, t1, jpuff2, t2 in zip(jpuffs_dt1[:-1], ts_dt1[:-1], jpuffs_dt1[1:], ts_dt1[1:]):  # set  = [time, state, idx]
        ts_dt1_new.append(t1)
        jpuffs_dt1_new.append(jpuff1)
        ts_dt1_new.append(t2)
        jpuffs_dt1_new.append(jpuff1)
    mean = np.mean(jpuffs)

    ax_ts2.set_xlabel("$t$")
    ax_ts2.set_xlim([0, 50])
    ax_ts2.set_ylabel(r"$y_{1.0}$")
    ax_ts2.set_ylim([0, 6])

    ax_ts2.plot(ts_dt1_new[plt_start:plt_stop], jpuffs_dt1_new[plt_start:plt_stop], c="C0")
    ax_ts2.axhline(N*mean_x, c="k", label="$\mu = {:.2f}$".format(mean))

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
    ax_hist1.plot(gauss001, fc.gaussian_dist(gauss001, mean, std01), c="C3")
    ax_hist1.hist(jpuffs_dt01, bins=50, alpha=0.7, color="C3", density=True,
                  label=r"$\Delta t = {:.1f}$".format(0.1) + "\n" + "$\sigma_{y_{0.1}}$" + " = {:.2f}".format(std01))
    legend_hs1 = ax_hist1.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_hs1.get_frame().set_linewidth(0.5)


    ax_hist2.set_xlabel("$y_{1.0}$")
    ax_hist2.set_xlim([-1, 6])
    ax_hist2.set_ylabel("$P(y_{1.0})$")
    ax_hist2.set_ylim([0, 1.25])
    ax_hist2.plot(gauss01, fc.gaussian_dist(gauss01, mean, std1), c="C0")
    ax_hist2.hist(jpuffs_dt1, bins=50, alpha=0.7, color="C0", density=True,
                  label=r"$\Delta t = {:.1f}$".format(1.0) + "\n" + "$\sigma_{y_{1.0}}$" + " = {:.2f}".format(std1))
    legend_hs2 = ax_hist2.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_hs2.get_frame().set_linewidth(0.5)

    plt.savefig(home + "/Data/Calcium/Plots/puff_current_coarse_grain.pdf", transparent=True)
    plt.show()
