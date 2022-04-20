import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

import styles as st


def get_mean_opn_channel(data):
    # data = [time, state, idx]
    data_single_cluster = []
    for set in data:
        if set[2] == 0:
            data_single_cluster.append(set)
    data_single_cluster_no_ref = []
    for set in data_single_cluster:
        if set[1] < 0:
            data_single_cluster_no_ref.append([set[0], 0, set[2]])
        else:
            data_single_cluster_no_ref.append(set)
    A = 0
    for set1, set2 in zip(data_single_cluster_no_ref[:-1], data_single_cluster_no_ref[1:]):
        time = set2[0] - set1[0]
        state = set1[1]
        A += time*state
    return A/1000


def get_var_open_channel(data):
    # data = [time, state, idx]
    data_single_cluster = []
    for set in data:
        if set[2] == 0:
            data_single_cluster.append(set)
    data_single_cluster_no_ref = []
    for set in data_single_cluster:
        if set[1] < 0:
            data_single_cluster_no_ref.append([set[0], 0, set[2]])
        else:
            data_single_cluster_no_ref.append(set)
    A2 = 0
    for set1, set2 in zip(data_single_cluster_no_ref[:-1], data_single_cluster_no_ref[1:]):
        time = set2[0] - set1[0]
        state = set1[1]
        A2 += time*(state**2)
    return A2/1000


def r_opn_full(x, y, n, a, b, c, K0, rmax):
    K = K0*(y**3 / (1 + y**3))
    r_opn_full = n * rmax * (np.power(x, a) / (1. + np.power(x, a))) * (np.power(K, c) / (np.power(K, c) + np.power(x, c)))
    return r_opn_full


def r_opn(x, y, n, a, b, rmax):
    r_opn = n * rmax * (np.power(x, a)/(1. + np.power(x, a)))*(np.power(y, b)/(1. + np.power(y, b)))
    return r_opn


def r_ref(x, y, n, a, b, rmax):
    r_ref = n * rmax * (np.power(x, a) / (1. + np.power(x, a))) * (np.power(y, b) / (1. + np.power(y, b)))
    return r_ref


if __name__ == "__main__":
    home = os.path.expanduser("~")
    n = 5
    m = 4
    a = 3
    b = 3
    c = 3
    K = 10
    r0_opn = 0.13
    r0_ref = 1.30
    r_opn_max = r0_opn*((1. + np.power(0.33, 3))/np.power(0.33, 3))*(2/1)
    r_ref_max = r0_ref*((1. + np.power(0.33, 3))/np.power(0.33, 3))*(2/1)

    cas_sim = np.linspace(0.33 , 0.99, 67)
    cas_theo = np.linspace(0.1, 10, 300)

    mean_x0s_sim = []
    var_x0s_sim = []

    mean_x0s_sim_ns = []
    var_x0s_sim_ns = []
    mean_x0s_ns = []
    var_x0s_ns = []

    tau_opns = []
    tau_clss = []
    p_opn1s = []
    p_opn2s = []
    p_opn3s = []

    tau_opns_full = []
    tau_clss_full = []
    p_opn1s_full = []
    p_opn2s_full = []
    p_opn3s_full = []

    ip3s = [0.1, 1.0, 10.0]
    num_chas = [1, 3, 5]

    for num_cha in num_chas:
        num_clu = int(100/num_cha)
        mean_x0s_n = []
        var_x0s_n = []
        mean_x0s_sim_n = []
        var_x0s_sim_n = []
        for ca in cas_sim:
            print(ca)
            r_ref = num_cha * r_ref_max * (np.power(ca, 3) / (1. + np.power(ca, 3))) * (
                        np.power(1., 3) / (1. + np.power(1., 3)))
            r_opn = num_cha * r_opn_max * (np.power(ca, 3) / (1. + np.power(ca, 3))) * (
                        np.power(1., 3) / (1. + np.power(1., 3)))
            r_cls = 50.
            tau_opn = (num_cha + 1) / (2 * r_cls)
            tau_cls = 1. / r_opn + (m - 1) / (r_ref)
            p_opn = tau_opn / (tau_opn + tau_cls)
            mean_x0 = ((num_cha + 2) / 3) * p_opn
            mean_x0s_n.append(mean_x0)
            var_x0 = (p_opn * (num_cha + 1) * (num_cha + 2) / 6 - mean_x0 ** 2)
            var_x0s_n.append(var_x0)

            folder = home + f"/Data/calcium_spikes_markov/ca_fix/clustersize_n{num_cha:d}/"
            file = f"puff_markov_cafix{ca:.2f}_ip1.00_tau1.00e+00_j1.00e+00_N{num_clu:d}_{num_cha:d}.dat"
            data = np.loadtxt(folder + file)
            mean_x0_sim = get_mean_opn_channel(data)
            var_x0_sim = get_var_open_channel(data) - mean_x0_sim ** 2
            mean_x0s_sim_n.append(mean_x0_sim)
            var_x0s_sim_n.append(var_x0_sim)

        mean_x0s_ns.append(mean_x0s_n)
        var_x0s_ns.append(var_x0s_n)

        mean_x0s_sim_ns.append(mean_x0s_sim_n)
        var_x0s_sim_ns.append(var_x0s_sim_n)


    for ca in cas_theo:
        ip3 = 1.
        for i, ip3 in enumerate(ip3s):
            r_ref = n*r_ref_max*(np.power(ca, 3)/(1. + np.power(ca, 3)))*(np.power(ip3, 3)/(1. + np.power(ip3, 3)))
            r_opn = n*r_opn_max*(np.power(ca, 3)/(1. + np.power(ca, 3)))*(np.power(ip3, 3)/(1. + np.power(ip3, 3)))
            r_ref_full = n*r_ref_max*(np.power(ca, 3)/(1. + np.power(ca, 3)))*(np.power(ip3, 3)/(1. + np.power(ip3, 3)))*(5./(5. + np.power(ca, 3)))
            r_opn_full = n*r_opn_max*(np.power(ca, 3)/(1. + np.power(ca, 3)))*(np.power(ip3, 3)/(1. + np.power(ip3, 3)))*(5./(5. + np.power(ca, 3)))
            r_cls = 50.
            if i == 0:
                tau_opn = (n+1)/(2*r_cls)
                tau_cls = 1./r_opn + (m-1)/(r_ref)
                p_opn1s.append(tau_opn / (tau_opn + tau_cls))

                tau_opn_full = (n+1)/(2*r_cls)
                tau_cls_full = 1./r_opn_full + (m-1)/(r_ref_full)
                p_opn1s_full.append(tau_opn_full / (tau_opn_full + tau_cls_full))
            if i == 1:
                tau_opn = (n+1)/(2*r_cls)
                tau_cls = 1./r_opn + (m-1)/(r_ref)
                p_opn2 = tau_opn / (tau_opn + tau_cls)
                p_opn2s.append(p_opn2)
                tau_opns.append(tau_opn)
                tau_clss.append(tau_cls)

                tau_opn_full = (n+1)/(2*r_cls)
                tau_cls_full = 1./r_opn_full + (m-1)/(r_ref_full)
                p_opn2s_full.append(tau_opn_full / (tau_opn_full + tau_cls_full))
                tau_opns_full.append(tau_opn_full)
                tau_clss_full.append(tau_cls_full)
            if i == 2:
                tau_opn = (n + 1) / (2 * r_cls)
                tau_cls = 1. / r_opn + (m - 1) / (r_ref)
                p_opn3s.append(tau_opn / (tau_opn + tau_cls))

                tau_opn_full = (n+1)/(2*r_cls)
                tau_cls_full = 1./r_opn_full + (m-1)/(r_ref_full)
                p_opn3s_full.append(tau_opn_full / (tau_opn_full + tau_cls_full))

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(9, 3.0))
    gs = gridspec.GridSpec(nrows=2, ncols=4, width_ratios=[1, 1, 1, 1])
    ax1 = fig.add_subplot(gs[:, 0:2])
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[1, 2])
    ax4 = fig.add_subplot(gs[0, 3])
    ax5 = fig.add_subplot(gs[1, 3])
    #ax3 = fig.add_subplot(grids[2])
    axis = [ax1, ax2, ax3, ax4, ax5]
    st.remove_top_right_axis(axis)

    ax1.text(0.10, 0.95, "A$_{i}$", fontsize=15, transform=ax1.transAxes, va='top')
    ax2.text(0.15, 0.95, "A$_{ii}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.15, 0.95, "A$_{iii}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.15, 0.95, "B$_{i}$", fontsize=13, transform=ax4.transAxes, va='top')
    ax5.text(0.15, 0.95, "B$_{ii}$", fontsize=13, transform=ax5.transAxes, va='top')

    ax1.text(1.0, 0.80, "Firing threshold", bbox=dict(facecolor='w', edgecolor="C7", pad=5.0), fontsize=11, va="center", ha="center")
    ax1.text(0.8, 0.72, "$s = 10$", fontsize=11, transform=ax1.transAxes, va='top')
    ax1.text(0.8, 0.55, "$s = 1.0$", fontsize=11, transform=ax1.transAxes, va='top')
    ax1.text(0.55, 0.08, "$s = 0.1$", fontsize=11, transform=ax1.transAxes, va='top')

    ax1.set_xlabel("$c_i$")
    ax1.set_xticks([0., 0.5, 1.0, 1.5, 2.0])
    ax1.set_xscale("log")
    ax1.set_ylabel(r"$p_{\rm open}$")
    ax1.set_ylim([0, 1])
    ax1.set_xticks([0.1, 0.33, 1, 10])
    ax1.set_xticklabels(["$10^{-1}$", "$c_R$", "$c_T = 10^0$", "$10^1$"])
    ax1.axvline(1, color="C7", alpha=0.7, zorder=1)
    ax1.plot(cas_theo, p_opn1s, color=st.colors[0])
    ax1.plot(cas_theo, p_opn2s, color=st.colors[0])
    ax1.plot(cas_theo, p_opn3s, color=st.colors[0])
    ax1.plot(cas_theo, p_opn1s_full, color=st.colors[0], ls=":")
    ax1.plot(cas_theo, p_opn2s_full, color=st.colors[0], ls=":")
    ax1.plot(cas_theo, p_opn3s_full, color=st.colors[0], ls=":")
    ax1.axvspan(1, 10, alpha=0.3, color="C7")

    ax2.set_xlabel("$c_i$")
    ax2.set_ylabel(r"$\tau_{\rm open}$  / s")
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.axvline(1, color="C7", alpha=0.7, zorder=1)
    ax2.axvspan(1, 10, alpha=0.3, color="C7")
    ax2.plot(cas_theo, tau_opns, color=st.colors[0])
    ax2.plot(cas_theo, tau_opns_full, color=st.colors[0], ls=":")

    ax3.set_xlabel("$c_i$")
    ax3.set_ylabel(r"$\tau_{\rm closed}$ / s")
    ax3.set_yscale("log")
    ax3.set_xscale("log")
    ax3.axvline(1, color="C7", alpha=0.7, zorder=1)
    ax3.axvspan(1, 10, alpha=0.3, color="C7")
    ax3.plot(cas_theo, tau_clss, color=st.colors[0])
    ax3.plot(cas_theo, tau_clss_full, color=st.colors[0], ls=":")

    ax4.set_xlim([0.33, 1.25])
    ax4.set_xticks([0.33, 1])
    ax4.set_ylim([0, 1.2*20*mean_x0s_ns[2][-1]])
    ax4.axvline(1, color="C7", alpha=0.7, zorder=1)
    ax4.axvspan(1, 10, alpha=0.3, color="C7")
    ax4.set_xticklabels(["$c_R$", "$c_T$"])
    ax4.set_xlabel("$c_i$")
    ax4.set_ylabel(r"$\langle j_{\rm puff} \rangle$")
    ax4.plot(cas_sim, [100*mean for mean in mean_x0s_sim_ns[0]], color=st.colors[0])
    ax4.plot(cas_sim, [33*mean for mean in mean_x0s_sim_ns[1]], color=st.colors[0])
    ax4.plot(cas_sim, [20*mean for mean in mean_x0s_sim_ns[2]], color=st.colors[0])
    #ax4.plot(cas_sim, mean_x0s, color=st.colors[2])
    ax4.plot(cas_sim, [100*mean for mean in mean_x0s_ns[0]], color=st.colors[2])
    ax4.plot(cas_sim, [33*mean for mean in mean_x0s_ns[1]], color=st.colors[2])
    ax4.plot(cas_sim, [20*mean for mean in mean_x0s_ns[2]], color=st.colors[2])
    ax4.text(1.01, 20*mean_x0s_ns[2][-1], "$n = 5$", fontsize=9, va='center')
    ax4.text(1.01, 33*mean_x0s_ns[1][-1], "$n = 3$", fontsize=9, va='center')
    ax4.text(1.01, 100*mean_x0s_ns[0][-1], "$n = 1$", fontsize=9, va='center')

    ax5.set_xlim([0.33, 1.25])
    ax5.set_xticks([0.33, 1])
    ax5.set_ylim([0, 1.2*np.sqrt(20*var_x0s_ns[2][-1])])
    ax5.axvline(1, color="C7", alpha=0.7, zorder=1)
    ax5.axvspan(1, 10, alpha=0.3, color="C7")
    ax5.set_xticklabels(["$c_R$", "$c_T$"])
    ax5.set_xlabel("$c_i$")
    ax5.set_ylabel(r"$\sqrt{\langle \Delta j_{\rm puff}^2 \rangle}$")
    ax5.plot(cas_sim, [np.sqrt(100*var) for var in var_x0s_sim_ns[0]], color=st.colors[0])
    ax5.plot(cas_sim, [np.sqrt(33*var) for var in var_x0s_sim_ns[1]], color=st.colors[0])
    ax5.plot(cas_sim, [np.sqrt(20*var) for var in var_x0s_sim_ns[2]], color=st.colors[0])
    #ax5.plot(cas_sim, [np.sqrt(var)/mean for var, mean in zip(var_x0s, mean_x0s)], color=st.colors[2])
    ax5.plot(cas_sim, [np.sqrt(100*var) for var in var_x0s_ns[0]], color=st.colors[2])
    ax5.plot(cas_sim, [np.sqrt(33*var) for var in var_x0s_ns[1]], color=st.colors[2])
    ax5.plot(cas_sim, [np.sqrt(20*var) for var  in var_x0s_ns[2]], color=st.colors[2])

    ax5.text(1.01, np.sqrt(20*var_x0s_ns[2][-1]), "$n = 5$", fontsize=9, va='center')
    ax5.text(1.01, np.sqrt(33*var_x0s_ns[1][-1]), "$n = 3$", fontsize=9, va='center')
    ax5.text(1.01, np.sqrt(100*var_x0s_ns[0][-1]), "$n = 1$", fontsize=9, va='center')
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig3.pdf",transparent=True)
    plt.show()