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

def get_As_by_n0(data, n):
    As = [[] for _ in range(n)]
    A_puff = 0
    puffing = False
    for ele1, ele2 in zip(data[:-1], data[1:]):
        if ele1[1] > 0:
            if puffing == False:
                n0 = int(ele1[1])
                puffing = True
            A_puff += ele1[1] * (ele2[0] - ele1[0])
        if ele1[1] <= 0:
            puffing = False
            if A_puff != 0:
                As[n0-1].append(A_puff)
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
    fig = plt.figure(tight_layout=True, figsize=(9, 4.5))
    gs = gridspec.GridSpec(nrows=3, ncols=4, width_ratios=[1,1,1,1])
    #gs0 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    #gs1 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0])
    #gs2 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1])
    ax_a1 = fig.add_subplot(gs[0, 0:2])
    ax_a2 = fig.add_subplot(gs[1, 0])
    ax_a3 = fig.add_subplot(gs[1, 1])
    ax_a4 = fig.add_subplot(gs[2, 0])
    ax_a5 = fig.add_subplot(gs[2, 1])
    ax_b1 = fig.add_subplot(gs[0, 2:4])
    ax_b2 = fig.add_subplot(gs[1, 2])
    ax_b3 =  fig.add_subplot(gs[1, 3])
    ax_b4 = fig.add_subplot(gs[2, 2])
    ax_b5 =  fig.add_subplot(gs[2, 3])
    axes = [ax_a1, ax_a2, ax_a3, ax_a4, ax_a5, ax_b1, ax_b2, ax_b3, ax_b4, ax_b5]
    st.remove_top_right_axis(axes)
    home = os.path.expanduser("~")

    ax_a1.text(0.1, 0.95, "A$_{i}$", fontsize=13, transform=ax_a1.transAxes, va='top')
    ax_a2.text(0.1, 0.95, "A$_{ii}$", fontsize=13, transform=ax_a2.transAxes, va='top')
    ax_a3.text(0.1, 0.95, "A$_{iii}$", fontsize=13, transform=ax_a3.transAxes, va='top')
    ax_a4.text(0.1, 0.95, "A$_{iv}$", fontsize=13, transform=ax_a4.transAxes, va='top')
    ax_a5.text(0.1, 0.95, "A$_{v}$", fontsize=13, transform=ax_a5.transAxes, va='top')

    ax_b1.text(0.085, 0.95, "B$_{i}$", fontsize=13, transform=ax_b1.transAxes, va='top')
    ax_b2.text(0.2, 0.95, "B$_{ii}$", fontsize=13, transform=ax_b2.transAxes, va='top')
    ax_b3.text(0.1, 0.95, "B$_{iii}$", fontsize=13, transform=ax_b3.transAxes, va='top')
    ax_b4.text(0.2, 0.95, "B$_{iv}$", fontsize=13, transform=ax_b4.transAxes, va='top')
    ax_b5.text(0.1, 0.95, "B$_{v}$", fontsize=13, transform=ax_b5.transAxes, va='top')
    #ax1.set_title("A) Puff statistics", fontsize=14, loc="left")
    #ax5.set_title("B) Interpuff interval statistics", fontsize=13, loc="left")

    # Subplot 1: x(t) over t
    ax_a1.set_xlabel("$t$ / s")
    ax_a1.set_ylabel("$x_i(t)$")
    folder = "/Data/calcium_spikes_markov/ca_fix/old/"
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

    t_left = starts_ipi[3] - ts[0]
    t_right = stops_ipi[5] - ts[0]

    ax_a1.plot(ts[:500] - ts[0], xs[:500], lw=1, color=st.colors[0])
    ax_a1.set_xlim([t_left - 0.15, t_right + 0.15])
    ax_a1.fill_between(ts[:500] - ts[0], xs[:500], color=st.colors[0], alpha=0.5)

    ax_a2.set_ylim([0, 1 / 5 + 0.15])
    ax_a2.set_xlabel(r"$n_{0}$")
    ax_a2.set_ylabel(r"$P(n_{0})$")
    ns_open = get_ns_open(data)
    ns_over_i = [0] * 5
    for n in ns_open:
        ns_over_i[n - 1] += 1
    for i, n in enumerate(ns_over_i):
        ax_a2.bar(i + 1, n / len(ns_open), 0.8, align="center", color=st.colors[0])
    ax_a2.plot([0.6, 0.6, 5.4, 5.4], [0, 1 / 5, 1 / 5, 0], color=st.colors[2])
    ax_a2.set_xticks([1, 2, 3, 4, 5])
    ax_a2.set_yticks([0, 0.1, 0.2, 0.3])
    # ax1.plot([0, 1, 1, 2, 3, 4, 5, 5, 6], [0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0], lw=1, color=st.colors[2])

    ax_a3.set_xlabel("$A$")
    ax_a3.set_ylabel("$P(A)$")
    As = get_As(data_no_ref)
    ax_a3.hist(As, bins=25, color=st.colors[0], density=True)


    Ass = []
    Ass_by_n = []
    num_chas = np.arange(1, 7)
    for num_cha in num_chas:
        file = f"puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N1_{num_cha:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        As = get_As(data_n)
        Ass.append(As)
        As_by_n = get_As_by_n0(data_n, num_cha)
        Ass_by_n.append(As_by_n)

    for n in range(1, 7):
        for n0 in range(1, n):
            print(f"mean {n0:d}:", np.mean(Ass_by_n[n-1][n0-1]), n0*(n0+1)/100)
            print("var:", np.var(Ass_by_n[n-1][n0-1]), n0*(n0+1)*(2*n0+1)/(6*50**2))
        print(np.mean(Ass[n-1]), (n+1)*(n+2)/(6*50))
        print(np.var(Ass[n-1]), (n+1)*(n+1)*(n+2)*(n+2)/((6*50)**2))

    ax_a4.set_xlabel("$N$")
    ax_a4.set_ylabel(r"$\langle A \rangle$")
    ax_a4.scatter(num_chas, [np.mean(As) for As in Ass], ec=st.colors[0], fc="w", s=15, zorder=2)
    Ass_theo = [(1 / 50) * (n + 1) * (n + 2) / 6 for n in num_chas]
    ax_a4.plot(num_chas, Ass_theo, c=st.colors[2], zorder=1)
    ax_a4.set_xticks([1, 2, 3, 4, 5, 6])
    ax_a4.set_yticks([0, 0.1, 0.2])

    ax_a5.set_xlabel("$N$")
    ax_a5.set_ylabel(r"$\sqrt{\langle \Delta A^2 \rangle}$")
    ax_a5.scatter(num_chas, [np.std(As) for As in Ass], ec=st.colors[0], fc="w", s=15, zorder=2)
    #Ass_theo = [(1 / 50) * (n + 1) * (n + 2) / 6 for n in num_chas]
    #ax_a4.plot(num_chas, Ass_theo, c=st.colors[2], zorder=1)
    ax_a5.set_xticks([1, 2, 3, 4, 5, 6])
    ax_a5.set_yticks([0.0, 0.1, 0.2])

    # Subplot 1: x(t) over t
    ax_b1.set_xlabel("$t$ / s")
    ax_b1.set_ylabel("$x_i(t)$")
    folder = "/Data/calcium_spikes_markov/ca_fix/old/"
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

    t_left = starts_ipi[5] - ts[0]
    t_right = stops_ipi[5] - ts[0]
    dt = t_right - t_left
    ax_b1.arrow(t_left, 2.8, dt, 0, fc="k", length_includes_head=True, head_width=0.25, head_length=0.5,
                lw=1.0, clip_on=False)
    ax_b1.arrow(t_right, 2.8, -dt, 0, fc="k", length_includes_head=True, head_width=0.25,
                head_length=0.5, lw=1.0, clip_on=False)
    ax_b1.text((t_left + t_right) / 2, 3.5, "$I_i$", va="center", ha="center")

    ax_b1.plot(ts[:500] - ts[0], xs[:500], lw=1, color=st.colors[0])
    ax_b1.set_xlim([0, 25])


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

    ax_b2.set_xlabel("$I$ / s")
    ax_b2.set_ylabel("$P(I)$")
    ax_b2.hist(ipis, bins=50, color=st.colors[0], density=True)
    ax_b2.plot(ts, p_ipi, color=st.colors[2])
    ax_b2.set_xlim([0, 10])

    ax_b3.set_xlabel(r"$\langle I \rangle$ / s")
    ax_b3.set_ylabel(r"$\sqrt{\langle \Delta I^2 \rangle}$ / s")
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
        ax_b3.plot(xs, [CV * x for x in xs], c=st.colors[2], zorder=1)
        ax_b3.plot(xs, xs, lw=1, ls=":", c="k", zorder=1)
        for i in range(5):
            if 100*(i+1) > len(ipis):
                break
            else:
                ipis_tmp = ipis[100*i:100*(i+1)]
                ax_b3.scatter(np.mean(ipis_tmp), np.std(ipis_tmp), ec=st.colors[0], fc="w", s=15, zorder=2)
    ax_b3.set_xlim([0, 6])
    ax_b3.set_ylim([0, 6])
    style = "Simple, head_width=3, head_length=5"
    kw = dict(arrowstyle=style, color="k")
    arrow = patches.FancyArrowPatch((4.2, 4.2), (5, 0),
                             connectionstyle="arc3,rad=-0.25", zorder=5, **kw)
    ax_b3.add_patch(arrow)
    ax_b3.text(5.2, 1.5, "$M$")
    ts = np.linspace(0, 1, 101)

    Iss = []
    std_Is = []
    mean_Is = []
    CV_Is = []
    for num_cha in num_chas:
        file = f"puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N1_{num_cha:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        Is = get_ipis(data_n)
        Iss.append(Is)
        ropn = num_cha * (1. + ratio * (num_cls - 1.)) / mean_ibi
        rref = num_cha * (1. + ratio * (num_cls - 1.)) / (ratio * mean_ibi)
        mean = (num_cls - 1) / rref + 1 / ropn
        std = np.sqrt((num_cls - 1) * (1 / rref) ** 2 + (1 / ropn) ** 2)
        mean_Is.append(mean)
        std_Is.append(std)
        CV_Is.append(std/mean)

    ax_b4.set_xlabel("$N$")
    ax_b4.set_ylabel(r"$\langle I \rangle$")
    ax_b4.set_xticks([1, 2, 3, 4, 5, 6])
    ax_b4.scatter(num_chas, [np.mean(Is) for Is in Iss], ec=st.colors[0], fc="w", s=15, zorder=2)
    ax_b4.plot(num_chas, mean_Is, c=st.colors[2], zorder=1)
    ax_b4.set_yticks([0, 5, 10])
    ax_b4.set_ylim([0, 11])

    ax_b5.set_xlabel("$N$")
    ax_b5.set_ylabel(r"$\sqrt{\langle \Delta I^2 \rangle}$")
    ax_b5.set_xticks([1, 2, 3, 4, 5, 6])
    ax_b5.scatter(num_chas, [np.std(Is) for Is in Iss], ec=st.colors[0], fc="w", s=15, zorder=2)
    ax_b5.plot(num_chas, std_Is, c=st.colors[2], zorder=1)
    ax_b5.set_yticks([0, 5, 10])
    ax_b5.set_ylim([0, 11])
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig2b.pdf", transparent=True)
    plt.show()
