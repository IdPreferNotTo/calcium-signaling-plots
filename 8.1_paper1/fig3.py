import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches
import os

import styles as st

def get_ns_open(data):
    ns_opn = []
    for ele1, ele2, in zip(data[:-1], data[1:]):
        if ele1[1] == 0:
           ns_opn.append(int(ele2[1]))
    return ns_opn


def get_As(data):
    # data = [t, state, idx]
    # Turn data into release ca per puff
    As = []
    A_puff = 0
    for ele1, ele2 in zip(data[:-1], data[1:]):
        if ele1[1] > 0:
            A_puff += ele1[1]*(ele2[0] - ele1[0])
        if ele1[1] <= 0:
            if A_puff != 0:
                As.append(A_puff)
                A_puff = 0
    return As


def get_ipis(data):
    is_ipi = False
    t0 = data[0][0]
    ipis= []
    for elem in data:
        if elem[1] < 0 and not is_ipi:
            is_ipi = True
            t0 = elem[0]
        if elem[1] > 0 and is_ipi:
            ipi = elem[0] - t0
            ipis.append(ipi)
            is_ipi = False
    return ipis


if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(9, 3))
    gs = gridspec.GridSpec(nrows=2, ncols=4, width_ratios=[1,1,1,1])
    #gs0 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    #gs1 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0])
    #gs2 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[0, 2:4])
    ax6 = fig.add_subplot(gs[1, 2])
    ax7 =  fig.add_subplot(gs[1, 3])
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
    st.remove_top_right_axis(axes)
    home = os.path.expanduser("~")

    ax1.text(0.1, 0.95, "A$_{i}$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, "A$_{ii}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.1, 0.95, "A$_{iii}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.1, 0.95, "A$_{iv}$", fontsize=13, transform=ax4.transAxes, va='top')

    ax5.text(0.085, 0.95, "B$_{i}$", fontsize=13, transform=ax5.transAxes, va='top')
    ax6.text(0.2, 0.95, "B$_{ii}$", fontsize=13, transform=ax6.transAxes, va='top')
    ax7.text(0.1, 0.95, "B$_{iii}$", fontsize=13, transform=ax7.transAxes, va='top')
    #ax1.set_title("A) Puff statistics", fontsize=14, loc="left")
    #ax5.set_title("B) Interpuff interval statistics", fontsize=13, loc="left")

    # Subplot 1: x(t) over t
    ax1.set_xlabel("$t$ / s")
    ax1.set_ylabel("$x_i(t)$")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/ca_fix/"
    data = np.loadtxt(home + folder + "puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N1_5.dat")
    data_tmp = []

    data_no_ref = [(t, x, i) if x >= 0 else (t, 0, i) for (t, x, i) in data]
    data2 = []
    for set1, set2 in zip(data_no_ref[:-1], data_no_ref[1:]):
        data2.append(set1)
        data2.append([set2[0], set1[1], set1[2]])
    ts, xs, idxs = np.transpose(data2)

    starts_ipi = []
    stops_ipi = []
    is_puff = False
    for elem in data:
        if elem[1] > 0 and not is_puff:
            is_puff = True
            t_start = elem[0]
            starts_ipi.append(t_start)
        if elem[1] < 0 and is_puff:
            is_puff = False
            t_stop = elem[0]
            stops_ipi.append(t_stop)

    t_left = starts_ipi[8] - ts[0]
    t_right = stops_ipi[8] - ts[0]

    ax1.plot(ts[:500] - ts[0], xs[:500], lw=1, color=st.colors[0])
    ax1.set_xlim([t_left - 0.15, t_right + 0.15])
    ax1.fill_between(ts[:500] - ts[0], xs[:500], color=st.colors[0], alpha=0.5)

    ax2.set_ylim([0, 1 / 5 + 0.15])
    ax2.set_xlabel(r"$n_{\rm opn}$")
    ax2.set_ylabel(r"$P(n_{\rm opn})$")
    ns_open = get_ns_open(data)
    ns_over_i = [0] * 5
    for n in ns_open:
        ns_over_i[n - 1] += 1
    for i, n in enumerate(ns_over_i):
        ax2.bar(i + 1, n / len(ns_open), 0.8, align="center", color=st.colors[0])
    ax2.set_xticks([1, 2, 3, 4, 5])
    ax2.set_yticks([0, 0.1, 0.2, 0.3])
    # ax1.plot([0, 1, 1, 2, 3, 4, 5, 5, 6], [0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0], lw=1, color=st.colors[2])

    ax3.set_xlabel("$A$")
    ax3.set_ylabel("$P(A)$")
    As = get_As(data_no_ref)
    ax3.hist(As, bins=25, color=st.colors[0], density=True)

    Ass = []
    num_chas = np.arange(1, 7)
    for num_cha in num_chas:
        file = f"puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N1_{num_cha:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        As = get_As(data_n)
        Ass.append(As)

    ax4.set_xlabel("$n$")
    ax4.set_ylabel(r"$\langle A \rangle$")
    ax4.scatter(num_chas, [np.mean(As) for As in Ass], ec=st.colors[0], fc="w", s=15, zorder=2)
    Ass_theo = [(1 / 50) * (n + 1) * (n + 2) / 6 for n in num_chas]
    ax4.plot(num_chas, Ass_theo, lw=1, c=st.colors[2], zorder=1)
    ax4.set_xticks([1, 2, 3, 4, 5, 6])

    # Subplot 1: x(t) over t
    ax5.set_xlabel("$t$ / s")
    ax5.set_ylabel("$x_i(t)$")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/ca_fix/"
    data = np.loadtxt(home + folder + "puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N1_5.dat")

    data_no_ref = [(t, x, i) if x>= 0 else (t, 0, i) for (t, x, i) in data]
    data_plot = []
    for set1, set2 in zip(data_no_ref[:-1], data_no_ref[1:]):
        data_plot.append(set1)
        data_plot.append([set2[0], set1[1], set1[2]])
    ts, xs, idxs = np.transpose(data_plot)

    ipis = get_ipis(data)
    starts_ipi =  []
    stops_ipi = []
    is_ipi = False
    for elem in data:
        if elem[1] < 0 and not is_ipi:
            is_ipi = True
            t_start = elem[0]
            starts_ipi.append(t_start)
        if elem[1] > 0 and is_ipi:
            is_ipi = False
            t_stop = elem[0]
            stops_ipi.append(t_stop)

    print(starts_ipi[8])
    t_left = starts_ipi[8] - ts[0]
    t_right = stops_ipi[8] - ts[0]
    dt = t_right - t_left
    ax5.arrow(t_left, 2.8, dt, 0, fc="k", length_includes_head=True, head_width=0.25, head_length=0.5,
              lw=1.0, clip_on=False)
    ax5.arrow(t_right, 2.8, -dt, 0, fc="k", length_includes_head=True, head_width=0.25,
              head_length=0.5, lw=1.0, clip_on=False)
    ax5.text((t_left + t_right) / 2, 3.5, "$I_i$", va="center", ha="center")

    ax5.plot(ts[:500] - ts[0], xs[:500], lw=1, color=st.colors[0])
    ax5.set_xlim([0, 25])


    # Subplot 2: IPI density
    num_cha = 5
    num_cls = 4
    mean_ibi = 10
    ratio = 0.1
    ropn = num_cha*(1. + ratio * (num_cls - 1.)) / mean_ibi
    rref = num_cha*(1. + ratio * (num_cls - 1.)) / (ratio * mean_ibi)
    mean = (num_cls-1)/rref + 1/ropn
    var = (num_cls-1)*(1/rref)**2 + (1/ropn)**2
    dr = rref - ropn
    ts = np.linspace(0, max(ipis), 100)
    p_ipi =[]
    for t in ts:
        c4 = (1 / 2) * np.power(rref, 3) * (np.exp(-dr * t) * (-dr * t * (dr * t + 2) - 2) + 2) / np.power(dr, 3)
        p4 = c4 * np.exp(-ropn * t)
        p_ipi.append(ropn * p4)

    print(np.mean(ipis), mean)
    ax6.set_xlabel("$I$ / s")
    ax6.set_ylabel("$P(I)$")
    ax6.hist(ipis, bins=50, color=st.colors[0], density=True)
    ax6.plot(ts, p_ipi, color=st.colors[2])
    ax6.set_xlim([0, 10])

    ax7.set_xlabel(r"$\langle I \rangle$ / s")
    ax7.set_ylabel(r"$\sqrt{\langle \Delta I^2 \rangle}$ / s")
    Ass = []
    num_chas = np.arange(1, 7)
    for num_cha in num_chas:
        file = f"puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N1_{num_cha:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        ipis = get_ipis(data_n)
        r = 0.1
        CV2 = (1. + (num_cls-1) * r**2)/(1. + (num_cls-1)*r)**2
        CV = np.sqrt(CV2)
        xs = np.linspace(0, 6, 100)
        ax7.plot(xs, [CV * x for x in xs], lw=1, c=st.colors[2], zorder=1)
        ax7.plot(xs, xs, lw=1, ls=":", c="k", zorder=1)
        for i in range(5):
            if 100*(i+1) > len(ipis):
                break
            else:
                ipis_tmp = ipis[100*i:100*(i+1)]
                ax7.scatter(np.mean(ipis_tmp), np.std(ipis_tmp), ec=st.colors[0], fc="w", s=15, zorder=2)
    ax7.set_xlim([0, 6])
    ax7.set_ylim([0, 6])
    style = "Simple, head_width=3, head_length=5"
    kw = dict(arrowstyle=style, color="k")
    arrow = patches.FancyArrowPatch((4.2, 4.2), (5, 0),
                             connectionstyle="arc3,rad=-0.25", zorder=5, **kw)
    ax7.add_patch(arrow)
    ax7.text(5.2, 1.5, "$m$")
    ts = np.linspace(0, 1, 101)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig3.png",
                transparent=True)
    plt.show()
