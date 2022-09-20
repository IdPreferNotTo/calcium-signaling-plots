import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, TransformedBbox, BboxPatch, BboxConnector
import os

import functions as fc
import styles as st


if __name__ == "__main__":
    # Prepare data
    home = os.path.expanduser("~")
    cas = np.linspace(0.01, 0.99, num=99)
    means_sim = []
    means_theory = []
    vars_sim = []
    vars_theory = []

    jpuffs_dt1_ca03 = []
    jpuffs_dt1_ca09 = []
    for ca in cas:
        print(f"{ca:.2f}")
        ca_fix = ca
        ca_res = 0.33
        N = 10
        n = 5
        m = 4
        jca = 1
        tau = 1
        folder = "/Data/calcium_spikes_markov/ca_fix/"
        file = "ca_markov_cafix{:.2f}_ip1.00_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca_fix)
        data = np.loadtxt(home + folder + file)
        ts, cas_sim, jpuffs, adaps = np.transpose(data)

        means_sim.append(np.mean(jpuffs))
        dt0 = 0.1
        f = 10
        cg_list = fc.coarse_grain_list(jpuffs, f)
        vars_sim.append(np.var(cg_list) * dt0 * f)

        r_opn = 0.13 * np.power(ca_fix / ca_res, 3) * ((1 + ca_res**3)/(1 + ca_fix**3))
        r_ref = 1.3 * np.power(ca_fix / ca_res, 3) * ((1 + ca_res**3)/(1 + ca_fix**3))
        r_cls = 50

        p0s = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)
        xs = fc.get_states(n, m)
        idxs = [i for i in range(n+m)]

        mean = sum([N*x*p for x, p in zip(xs, p0s)])
        means_theory.append(mean)

        D_theory = 0
        for k in idxs:
            sum_over_i = 0
            f_from_k_to = fc.f_from_k_invert_M(k, r_ref, r_opn, r_cls, n, m)
            for i in idxs:
                sum_over_i += xs[i]*f_from_k_to[i]
            D_theory += xs[k] * p0s[k]*sum_over_i
        vars_theory.append(2*N*D_theory)

        if ca == 0.3:
            jpuffs_dt1_ca03 = fc.coarse_grain_list(jpuffs, 10)
        if ca == 0.9:
            jpuffs_dt1_ca09 = fc.coarse_grain_list(jpuffs, 10)

    stds_sim = [np.sqrt(vars) for vars in vars_sim]
    stds_theory = [np.sqrt(vars) for vars in vars_theory]

    std_upper_theory = [mean + std for mean, std in zip(means_theory, stds_theory)]
    std_lower_theory = [mean - std for mean, std in zip(means_theory, stds_theory)]
    std_upper_sim= [mean + std for mean, std in zip(means_sim, stds_sim)]
    std_lower_sim = [mean - std for mean, std in zip(means_sim, stds_sim)]

    # Plot data
    st.set_default_plot_style()
    fig = plt.figure(tight_layout = True, figsize=(64/9, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    ax = fig.add_subplot(gs[0:2, 0:2])
    ax1 = fig.add_subplot(gs[0, 2])
    ax2 = fig.add_subplot(gs[1, 2])
    st.remove_top_right_axis([ax, ax1, ax2])

    # Main ax
    ax.set_xlabel(r"$c_i$")
    ax.set_xlim([0, 1])
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    ax.set_ylabel(r"$\mu, \sqrt{2D_N}$")
    ax.plot(cas, means_theory, c="k", ls="--", label=r"$\mu$ theory")
    ax.plot(cas, means_sim, c="k", label =r"$\mu$")
    ax.fill_between(cas, std_upper_sim, std_lower_sim, color="C0", alpha=0.55, label=r"$\sigma_{y_1}\sqrt{\Delta t}, \Delta t = 1$")
    ax.plot(cas, std_upper_theory, c="C0", ls="--", label=r"$\sqrt{2D_N}$ theory")
    ax.plot(cas, std_lower_theory, c="C0", ls="--")
    ax.axhline(0, lw="1", ls=":", c="k")
    legend = ax.legend(loc=2, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)

    # Histogram for ca = 0.9
    ax1.set_xlabel("$y_1$")
    ax1.set_xlim([0, 8])
    ax1.set_ylabel("$P(y_1)$")
    ax1.set_ylim([0, 2])
    ax1.hist(jpuffs_dt1_ca09, bins=50, color="C0", alpha=0.55, density=True,
                  label=r"[Ca\textsuperscript{2+}]$ = 0.9$")
    mean_ca09 = np.mean(jpuffs_dt1_ca09)
    std_ca09 = np.std(jpuffs_dt1_ca09)
    gauss_ca09 = np.linspace(mean_ca09 - 3 * std_ca09, mean_ca09 + 3 * std_ca09, 100)
    ax1.plot(gauss_ca09, fc.gaussian_dist(gauss_ca09, mean_ca09, std_ca09), c="C7")
    legend = ax1.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)

    # Histogram for ca = 0.3
    ax2.set_xlabel("$y_1$")
    ax2.set_xlim([0, 8])
    ax2.set_ylabel("$P(y_1)$")
    ax2.set_ylim([0, 2])
    ax2.hist(jpuffs_dt1_ca03, bins=50, color="C0", alpha=0.55, density=True,
                  label=r"[Ca\textsuperscript{2+}]$ = 0.3$")
    mean_ca03 = np.mean(jpuffs_dt1_ca03)
    std_ca03 = np.std(jpuffs_dt1_ca03)
    gauss_ca03 = np.linspace(mean_ca03 - 3 * std_ca03, mean_ca03 + 3 * std_ca03, 100)
    ax2.plot(gauss_ca03, fc.gaussian_dist(gauss_ca03, mean_ca03, std_ca03), c="C7")
    legend = ax2.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)

    # Save and show figure
    plt.savefig(home + "/Plots/puff_current_mean_var_static.pdf", transparent=True)
    plt.show()