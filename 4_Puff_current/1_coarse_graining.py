import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    # Set up Parameters
    ci = 0.20
    cT = 0.50
    cR = 0.20
    K = 10
    N = 5
    M = 4
    p = 1
    tau = 1

    # Load Data
    data = df.load_traces_fixed_ci_markov(tau, p, K=10, ca=ci)
    ts, cas, jps, adaps = np.transpose(data)

    r_opn = df.r_opn(ci, 1, N)
    r_ref = df.r_ref
    r_cls = df.r_cls

    # alternatively one could use an analytical expression for the mean puff current
    p_open = fc.p_open_cluster_theory(ci, N, M)
    mean_jp = p * K * fc.mean_puff_strength_cluster_theory(N) / fc.tau_total_cluster_theory(ci, N, M)

    jp_dt1 = jps
    jp_dt100 = fc.coarse_grain_list(jps, 100)

    # Set up plot style
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4.5, 4))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    plt_start = 0
    plt_stop = plt_start + 500

    ax1.set_xlabel("$t$ / s")
    ax1.set_xlim([0, 25])
    ax1.set_ylabel(r"$j_{\rm puff}$")
    ax1.set_ylim([0, 6])
    ax1.plot(ts, jps, c=st.colors[0])
    ax1.axhline(mean_jp, c=st.colors[2])

    ax3.plot([t for t in ts[::100]], jp_dt100, c="C0")
    ax3.axhline(mean_jp, c=st.colors[2])

    ax3.set_xlabel("$t$")
    ax3.set_xlim([0, 50])
    ax3.set_ylabel(r"$y_{1.0}$")
    ax3.set_ylim([0, 6])



    std01 = np.std(jp_dt1)
    std1 = np.std(jp_dt100)
    var01 = np.var(jp_dt1)
    var1 = np.var(jp_dt100)
    gauss001 = np.linspace(mean_jp - 3 * std01, mean_jp + 3 * std01, 100)
    gauss01 = np.linspace(mean_jp - 3 * std1, mean_jp + 3 * std1, 100)

    ax2.set_xlabel(r"$j_{\rm puff}$")
    ax2.set_xlim([-1, 6])
    ax2.set_ylabel(r"$p(j_{\rm puff})$")
    ax2.plot(gauss001, fc.gaussian_dist(gauss001, mean_jp, std01), c=st.colors[0])
    ax2.hist(jp_dt1, bins=50, alpha=0.7, color=st.colors[0], density=True,
             label=r"$\Delta t = {:.1f}$".format(0.1) + "\n" + "$\sigma_{y_{0.1}}$" + " = {:.2f}".format(std01))
    legend_hs1 = ax2.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_hs1.get_frame().set_linewidth(0.5)


    ax4.set_xlabel("$y_{1.0}$")
    ax4.set_xlim([-1, 6])
    ax4.set_ylabel("$P(y_{1.0})$")
    ax4.set_ylim([0, 1.25])
    ax4.plot(gauss01, fc.gaussian_dist(gauss01, mean_jp, std1), c="C0")
    ax4.hist(jp_dt100, bins=50, alpha=0.7, color="C0", density=True,
             label=r"$\Delta t = {:.1f}$".format(1.0) + "\n" + "$\sigma_{y_{1.0}}$" + " = {:.2f}".format(std1))
    legend_hs2 = ax4.legend(fancybox=False, edgecolor="k", framealpha=1.0)
    legend_hs2.get_frame().set_linewidth(0.5)

    plt.show()
