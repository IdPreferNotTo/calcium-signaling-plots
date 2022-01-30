import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    # Set up Parameters
    ca = 0.5
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

    Popen = np.power(1 + 2 * np.power(n * (n + 1), -1, dtype=float) * np.power(ca_rest / ca, 3) * (1 + ca ** 3) / (1 + ca_rest ** 3) * 10 / 0.02, -1)
    mean_jpuff = jca * (N * (n + 2) / 3) * Popen

    cgrain_factor1 = 1
    cgrain_factor2 = 10
    jpuffs_f1 = fc.coarse_grain_list(jpuffs, cgrain_factor1)
    jpuffs_f2 = fc.coarse_grain_list(jpuffs, cgrain_factor2)

    ts_f1 = [cgrain_factor1 * t for t in ts]
    ts_f2 = [cgrain_factor2 * t for t in ts]
    ts_f1_new = []
    ts_f2_new = []
    jpuffs_f1_new = []
    jpuffs_f2_new = []
    for jpuff1, t1, jpuff2, t2 in zip(jpuffs_f1[:-1], ts_f1[:-1], jpuffs_f1[1:], ts_f1[1:]):  # set  = [time, state, idx]
        ts_f1_new.append(t1)
        jpuffs_f1_new.append(jpuff1)
        ts_f1_new.append(t2)
        jpuffs_f1_new.append(jpuff1)

    for jpuff1, t1, jpuff2, t2 in zip(jpuffs_f2[:-1], ts_f2[:-1], jpuffs_f2[1:], ts_f2[1:]):  # set  = [time, state, idx]
        ts_f2_new.append(t1)
        jpuffs_f2_new.append(jpuff1)
        ts_f2_new.append(t2)
        jpuffs_f2_new.append(jpuff1)

    mean = np.mean(jpuffs)

    std_f1 = np.std(jpuffs_f1)
    std_f2 = np.std(jpuffs_f2)
    var_f1 = np.var(jpuffs_f1)
    var_f2 = np.var(jpuffs_f2)
    covar_f1 = fc.k_corr(jpuffs_f1, jpuffs_f1, k=1)
    covar_f2 = fc.k_corr(jpuffs_f2, jpuffs_f2, k=1)
    print(covar_f1, covar_f2)
    gauss001 = np.linspace(mean - 3 * std_f1, mean + 3 * std_f1, 100)
    gauss01 = np.linspace(mean - 3 * std_f2, mean + 3 * std_f2, 100)

    # Set up plot style
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 2))
    gs = gridspec.GridSpec(1, 2)
    gs1 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], wspace=0)
    gs2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0)
    ax0 = fig.add_subplot(gs1[0:2])
    ax1 = fig.add_subplot(gs1[2])
    ax2 = fig.add_subplot(gs2[0:2])
    ax3 = fig.add_subplot(gs2[2])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3])

    ylim1 = 5
    ylim2 = 5
    xlim1 = 10
    xlim2 = 100
    ax0.set_xlabel("$t$ / s")
    ax0.set_xlim([0, xlim1])
    ax0.set_ylabel(rf"$X(t; \Delta t = {0.1*cgrain_factor1:.1f})$")
    ax0.set_ylim([0, ylim1])
    ax0.plot(ts_f1_new, jpuffs_f1_new, c=st.colors[0])
    ax0.axhline(mean_jpuff, c=st.colors[2])

    ax2.set_xlabel("$t$ / s")
    ax2.set_xlim([0, xlim2])
    ax2.set_ylabel(rf"$X(t; \Delta t = {0.1*cgrain_factor2:.1f})$")
    ax2.set_ylim([0, ylim2])
    ax2.plot(ts_f2_new, jpuffs_f2_new, c=st.colors[0])
    ax2.axhline(mean_jpuff, c=st.colors[2])

    ax1.set_xlabel("$P(X)$")
    ax1.set_ylim([0, 6.])
    ax1.set_yticks([])
    ax1.set_xticks([])
    ax1.plot(fc.gaussian_dist(gauss001, mean, std_f1), gauss001, c="k", lw=1.)
    ax1.hist(jpuffs_f1, bins=25, alpha=0.7, color=st.colors[0], density=True, orientation="horizontal")

    ax3.set_xlabel("$P(X)$")
    ax3.set_ylim([0, 6.])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.plot(fc.gaussian_dist(gauss01, mean, std_f2), gauss01, c="k", lw=1.)
    ax3.hist(jpuffs_f2, bins=7, alpha=0.7, color=st.colors[0], density=True, orientation="horizontal")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig9.png", transparent=True)
    plt.show()
