import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

import styles as st
import functions as fc


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

if __name__ == "__main__":
    home = os.path.expanduser("~")
    n = 5
    m = 3
    a = 3
    b = 3
    c = 3
    K = 10
    r0_opn = 0.1
    r_ref = 20.0
    r_cls = 50.
    r_opn_max = r0_opn*((1. + np.power(0.2, 3))/np.power(0.2, 3))*(2/1)

    cas_theo = np.linspace(0.05, 5, 300)

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
    for ca in cas_theo:
        ip3 = 1.
        for i, ip3 in enumerate(ip3s):
            r_opn = n*r_opn_max*(np.power(ca, 3)/(1. + np.power(ca, 3)))*(np.power(ip3, 3)/(1. + np.power(ip3, 3)))
            r_opn_full = n*r_opn_max*(np.power(ca, 3)/(1. + np.power(ca, 3)))*(np.power(ip3, 3)/(1. + np.power(ip3, 3)))*(10./(10. + np.power(ca, 3)))
            r_cls = 50.
            if i == 0:
                tau_opn = (n+1)/(2*r_cls)
                tau_cls = 1./r_opn + (m-1)/(r_ref)
                p_opn1s.append(tau_opn / (tau_opn + tau_cls))

                tau_opn_full = (n+1)/(2*r_cls)
                tau_cls_full = 1./r_opn_full + (m-1)/(r_ref)
                p_opn1s_full.append(tau_opn_full / (tau_opn_full + tau_cls_full))
            if i == 1:
                tau_opn = (n+1)/(2*r_cls)
                tau_cls = 1./r_opn + (m-1)/(r_ref)
                p_opn2 = tau_opn / (tau_opn + tau_cls)
                p_opn2s.append(p_opn2)
                tau_opns.append(tau_opn)
                tau_clss.append(tau_cls)

                tau_opn_full = (n+1)/(2*r_cls)
                tau_cls_full = 1./r_opn_full + (m-1)/(r_ref)
                p_opn2s_full.append(tau_opn_full / (tau_opn_full + tau_cls_full))
                tau_opns_full.append(tau_opn_full)
                tau_clss_full.append(tau_cls_full)
            if i == 2:
                tau_opn = (n + 1) / (2 * r_cls)
                tau_cls = 1. / r_opn + (m - 1) / (r_ref)
                p_opn3s.append(tau_opn / (tau_opn + tau_cls))

                tau_opn_full = (n+1)/(2*r_cls)
                tau_cls_full = 1./r_opn_full + (m-1)/(r_ref)
                p_opn3s_full.append(tau_opn_full / (tau_opn_full + tau_cls_full))

    ca_exp1 = np.logspace(-1, 0, 10)
    ca_exp1 = [x/2 for x in ca_exp1]
    p_open_exps = []
    tau_closeds = []
    tau_opens = []
    folder =  "/Data/calcium_spikes_markov/ca_fix/"
    for ip3 in [0.1, 1.0, 10]:
        p_open_exp = []
        tau_closed = []
        tau_open = []
        for ca in ca_exp1:
            file = f"puff_markov_cafix{ca:.2f}_ip{ip3:.2f}_tau1.00e+00_j1.00e+00_K1_5.dat"
            data = np.loadtxt(home + folder + file)
            if len(data) == 0:
                p_open = 0
                t_cls = 1000
                t_opn = 0
            else:
                p_open = fc.p_open_cluster_data(data)
                t_cls = fc.tau_close_cluster_data(data)
                t_opn = fc.tau_open_cluster_data(data)
            p_open_exp.append(p_open)
            tau_closed.append(t_cls)
            tau_open.append(t_opn)
        p_open_exps.append(p_open_exp)
        tau_closeds.append(tau_closed)
        tau_opens.append(tau_open)
    
    colors = st.Colors()
    st.set_default_plot_style()
    w = 3.25*1.25
    h = 3.0*1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(nrows=3, ncols=2)
    ax1 = fig.add_subplot(gs[0:2, :])
    ax2 = fig.add_subplot(gs[2, 0])
    ax3 = fig.add_subplot(gs[2, 1])
    axis = [ax1, ax2, ax3]
    st.remove_top_right_axis(axis)

    ax1.text(0.10, 0.95, "A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.15, 0.95, "B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.15, 0.95, "C", fontsize=11, transform=ax3.transAxes, va='top')

    #ax1.text(1.0, 0.91, "Firing threshold", bbox=dict(facecolor='w', edgecolor="C7", pad=5.0), fontsize=11, va="top", ha="center")
    ax1.text(0.55, 0.80, "$s = 10$", fontsize=10, transform=ax1.transAxes, va='top')
    ax1.text(0.55, 0.45, "$s = 1.0$", fontsize=10, transform=ax1.transAxes, va='top')
    ax1.text(0.55, 0.08, "$s = 0.1$", fontsize=10, transform=ax1.transAxes, va='top')

    ax1.set_xlabel(r"$c_{\rm i}$")
    ax1.set_xscale("log")
    ax1.set_ylabel(r"$p_{\rm opn}$")
    ax1.set_ylim([0, 0.5])
    ax1.set_xticks([0.05, 0.2, 0.5, 5])
    ax1.set_xticklabels(["$5\cdot 10^{-2}$", "$c_R$", "$c_T = 5 \cdot 10^{-1}$", "$5 \cdot 10^0$"])
    ax1.axvline(0.5, color="C7", alpha=0.7, zorder=1)

    ax1.scatter(ca_exp1, p_open_exps[0], s=25, fc="w", ec=colors.palette[0], zorder=3, label="Sim.")
    ax1.scatter(ca_exp1, p_open_exps[1], s=25, fc="w", ec=colors.palette[0], zorder=3)
    ax1.scatter(ca_exp1, p_open_exps[2], s=25, fc="w", ec=colors.palette[0], zorder=3)
    ax1.plot(cas_theo, p_opn1s, color=colors.palette[6], label="Theory", zorder=2)
    ax1.plot(cas_theo, p_opn2s, color=colors.palette[6], zorder=2)
    ax1.plot(cas_theo, p_opn3s, color=colors.palette[6], zorder=2)
    ax1.plot(cas_theo, p_opn1s_full, color=colors.palette[6], ls=":", zorder=2)
    ax1.plot(cas_theo, p_opn2s_full, color=colors.palette[6], ls=":", zorder=2)
    ax1.plot(cas_theo, p_opn3s_full, color=colors.palette[6], ls=":", zorder=2)

    ax1.set_xlim([0.05, 5])
    ax1.axvspan(0.5, 5, alpha=0.3, color="C7")
    ax1.legend(fancybox=False, fontsize=9, loc=1)

    ax2.set_xlabel(r"$c_{\rm i}$")
    ax2.set_ylabel(r"$\tau_{\rm opn}$  / s")
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylim([0.01, 0.15])
    ax2.axvline(0.5, color="C7", alpha=0.7, zorder=1)
    ax2.axvspan(0.5, 5, alpha=0.3, color="C7")
    ax2.plot(cas_theo, tau_opns, color=colors.palette[6])
    ax2.plot(cas_theo, tau_opns_full, color=colors.palette[6], ls=":")
    ax2.scatter(ca_exp1[1:], tau_opens[1][1:], s=15, fc="w", ec=colors.palette[0], zorder=2)

    ax3.set_xlabel(r"$c_{\rm i}$")
    ax3.set_ylabel(r"$\tau_{\rm cls}$ / s")
    ax3.set_yscale("log")
    ax3.set_xscale("log")
    ax3.axvline(0.5, color="C7", alpha=0.7, zorder=1)
    ax3.axvspan(0.5, 5, alpha=0.3, color="C7")
    ax3.plot(cas_theo, tau_clss, color=colors.palette[6])
    ax3.plot(cas_theo, tau_clss_full, color=colors.palette[6], ls=":")
    ax3.scatter(ca_exp1[1:], tau_closeds[1][1:], s=15, fc="w", ec=colors.palette[0], zorder=2)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/SUB2/figures/fig3.pdf",transparent=True)
    plt.show()