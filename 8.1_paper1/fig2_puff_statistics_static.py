import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle
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
            A_puff += ele1[1] * (ele2[0] - ele1[0])
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
                As[n0 - 1].append(A_puff)
                A_puff = 0
    return As


def get_ipis(data):
    is_ipi = False
    t0 = data[0][0]
    ipis = []
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
    gs = gridspec.GridSpec(nrows=3, ncols=4, width_ratios=[1, 1, 1, 1])
    # gs0 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    # gs1 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0])
    # gs2 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1])
    ax_a1 = fig.add_subplot(gs[0, 0:2])
    ax_a2 = fig.add_subplot(gs[1, 0])
    ax_a3 = fig.add_subplot(gs[1, 1])
    ax_a4 = fig.add_subplot(gs[2, 0])
    ax_a5 = fig.add_subplot(gs[2, 1])
    ax_b1 = fig.add_subplot(gs[0, 2:4])
    ax_b2 = fig.add_subplot(gs[1, 2])
    ax_b3 = fig.add_subplot(gs[1, 3])
    ax_b4 = fig.add_subplot(gs[2, 2])
    ax_b5 = fig.add_subplot(gs[2, 3])
    axes = [ax_a1, ax_a2, ax_a3, ax_a4, ax_a5, ax_b1, ax_b2, ax_b3, ax_b4, ax_b5]
    st.remove_top_right_axis(axes)
    colors = st.Colors()
    home = os.path.expanduser("~")

    ax_a1.text(0.05, 1.00, r"A$_{\rm i}$", fontsize=13, transform=ax_a1.transAxes, va='top')
    ax_a2.text(0.15, 1.00, r"A$_{\rm ii}$", fontsize=13, transform=ax_a2.transAxes, va='top')
    ax_a3.text(0.15, 1.00, r"A$_{\rm iii}$", fontsize=13, transform=ax_a3.transAxes, va='top')
    ax_a4.text(0.15, 1.00, r"A$_{\rm iv}$", fontsize=13, transform=ax_a4.transAxes, va='top')
    ax_a5.text(0.15, 1.00, r"A$_{\rm v}$", fontsize=13, transform=ax_a5.transAxes, va='top')

    ax_b1.text(0.085, 1.00, r"B$_{\rm i}$", fontsize=13, transform=ax_b1.transAxes, va='top')
    ax_b2.text(0.15, 1.00, r"B$_{\rm ii}$", fontsize=13, transform=ax_b2.transAxes, va='top')
    ax_b3.text(0.15, 1.00, r"B$_{\rm iii}$", fontsize=13, transform=ax_b3.transAxes, va='top')
    ax_b4.text(0.15, 1.00, r"B$_{\rm iv}$", fontsize=13, transform=ax_b4.transAxes, va='top')
    ax_b5.text(0.15, 1.00, r"B$_{\rm v}$", fontsize=13, transform=ax_b5.transAxes, va='top')

    # Subplot 1: x(t) over t
    ax_a1.set_xlabel("$t$ / s")
    ax_a1.set_ylabel("$x(t)$")
    folder = "/Data/calcium_spikes_markov/ca_fix/"
    data = np.loadtxt(home + folder + "puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_K1_5.dat")
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
            t_ipi_start = elem[0]
            starts_ipi.append(t_ipi_start)
        if elem[1] < 0 and is_puff:
            is_puff = False
            t_ipi_stop = elem[0]
            stops_ipi.append(t_ipi_stop)

    t_start = 30
    t_stop = 60
    puff_index1 = 28
    puff_index2 = puff_index1 + 1
    t_left = starts_ipi[puff_index1] - t_start
    t_right = stops_ipi[puff_index2] - t_start
    for i, (t1, t2) in enumerate(zip(ts[:-1], ts[1:])):
        if t1 <= t_start and t2 > t_start:
            t_start_index = i
        if t1 <= t_stop and t2 > t_stop:
            t_stop_index = i+1
        if t1 <= (t_left + t_start) and t2 > (t_left + t_start):
            t_left_index = i
        if t1 <= (t_right + t_start) and t2 > (t_right + t_start):
            t_right_index = i+1

    print(t_start_index, ts[t_start_index])
    print(t_stop_index, ts[t_stop_index])
    print(t_left, t_left_index, ts[t_left_index])
    print(t_right, t_right_index, ts[t_right_index])

    Dt = t_right - t_left + 0.3
    # plt.rcParams['hatch.linewidth'] = 2
    ax_a1.axvspan(t_left - 0.05, t_right + 0.05, alpha=0.5, hatch='////', edgecolor=colors.green[1], facecolor="none")
    ax_b1.axvspan(t_left - 0.25, t_right + 0.25, alpha=0.5, hatch='/////', edgecolor=colors.green[1], facecolor="none")

    ax_a1.plot(ts[t_left_index:t_right_index] - t_start, xs[t_left_index:t_right_index], lw=1, color=colors.palette[0])
    ax_a1.fill_between(ts[t_left_index:t_right_index] - t_start, xs[t_left_index:t_right_index], color="#8B9FAF")
    ax_a1.set_xlim([t_left - 0.15, t_right + 0.15])
    ax_a1.set_ylim([0, 6])
    ax_a1.text(stops_ipi[puff_index1] - t_start + 0.02, 3.5, "$A_{i}$", va="center", ha="center")
    ax_a1.text(starts_ipi[puff_index2] - t_start + 0.05, 3.5, "$A_{i+1}$", va="center", ha="center")

    line1, = ax_a1.plot([-1, -1], color=colors.palette[0], label="Sim.")
    line2, = ax_a1.plot([-1, -1], color=colors.palette[5], label="Theory")
    legend = ax_a1.legend(fancybox=False, framealpha=1., edgecolor="k", ncol=2, handles=[line1, line2], fontsize=9,
                          loc=1, bbox_to_anchor=(1.0, 1.1))
    legend.get_frame().set_linewidth(0.5)

    ax_a2.set_ylim([0, 1 / 5 + 0.15])
    ax_a2.set_xlabel(r"$n_{0}$")
    ax_a2.set_ylabel(r"$p(n_{0} | N)$")
    ns_open = get_ns_open(data)
    ns_over_i = [0] * 5
    for n in ns_open:
        ns_over_i[n - 1] += 1
    for i, n in enumerate(ns_over_i):
        ax_a2.bar(i + 1, n / len(ns_open), 0.8, align="center", color=colors.palette[0])
    ax_a2.plot([0.6, 0.6, 5.4, 5.4], [0, 1 / 5, 1 / 5, 0], color=colors.palette[5])
    ax_a2.set_xticks([1, 2, 3, 4, 5])
    ax_a2.set_yticks([0, 0.1, 0.2, 0.3])
    # ax1.plot([0, 1, 1, 2, 3, 4, 5, 5, 6], [0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0], lw=1, color=colors.palette[5])

    ax_a3.set_xlabel("$A$")
    ax_a3.set_ylabel("$p(A)$")
    As = get_As(data_no_ref)
    ax_a3.hist(As, bins=25, color=colors.palette[0], density=True)

    Ass = []  # As are the Areas of puffs for a given number of total channels N.
    # Ass is the list that contains all those As for N=1 to N=6.
    Ass_over_n = []  # As_by_n are the Area of puffs given a number of responding channels n.
    # As_by_n is thus already a list of lists namely As_by_n = [(As | n=1), (As | n=2), ...]
    # Ass_by_n is the list that contains all those conditional Areas for varying numbers of total channels N.
    Ns = np.arange(1, 8)
    for N in Ns:
        file = f"puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_K1_{N:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        As = get_As(data_n)
        Ass.append(As)
        As_for_n = get_As_by_n0(data_n, N)
        Ass_over_n.append(As_for_n)

    for n in range(1, 7):
        for n0 in range(1, n + 1):
            print(f"mean {n0:d}:", np.mean(Ass_over_n[n - 1][n0 - 1]), n0 * (n0 + 1) / 100)
            print("var:", np.var(Ass_over_n[n - 1][n0 - 1]), n0 * (n0 + 1) * (2 * n0 + 1) / (6 * 50 ** 2))
        print(np.mean(Ass[n - 1]), (n + 1) * (n + 2) / (6 * 50))
        print(np.var(Ass[n - 1]), (n + 1) * (n + 1) * (n + 2) * (n + 2) / ((6 * 50) ** 2))
    Ass_theo = [(1 / 50) * (n + 1) * (n + 2) / 6 for n in Ns]

    Mean_of_VarA_over_n = []
    Mean_of_VarA_over_n_theo = []
    Var_of_MeanA_over_n = []
    Var_of_MeanA_over_n_theo = []

    rcls = 50.
    for n in range(1, 7):
        VarA_for_n = 0
        MeanA_over_n0 = []
        for n0 in range(1, n + 1):
            VarA_for_n += (1 / n) * np.power(6 * rcls * rcls, -1) * n0 * (n0 + 1) * (2 * n0 + 1)
            MeanA_for_n0 = np.power(2 * rcls, -1) * n0 * (n0 + 1)
            MeanA_over_n0.append(MeanA_for_n0)
        Var_of_MeanA_over_n.append(np.var(MeanA_over_n0))
        Var_of_MeanA_over_n_theo.append(np.mean([x ** 2 for x in MeanA_over_n0]) - (np.mean(MeanA_over_n0)) ** 2)
        Mean_of_VarA_over_n.append(VarA_for_n)
        Mean_of_VarA_over_n_theo.append(np.power(6 * rcls * rcls, -1) * (1 / 2) * (n + 1) * (n + 1) * (n + 2))
    print(Mean_of_VarA_over_n)
    print(Mean_of_VarA_over_n_theo)

    print(Var_of_MeanA_over_n)
    print(Var_of_MeanA_over_n_theo)

    ax_a4.set_xlabel("$N$")
    ax_a4.set_ylabel(r"$\langle A \rangle$")
    ax_a4.scatter(Ns, [np.mean(As) for As in Ass], ec=colors.palette[0], fc="w", s=15, zorder=2)
    ax_a4.plot(Ns, Ass_theo, c=colors.palette[5], zorder=1)
    ax_a4.set_xticks(Ns)
    ax_a4.set_yticks([0, 0.1, 0.2])

    ax_a5.set_xlabel("$N$")
    ax_a5.set_ylabel(r"${CV}_A^2$")
    ax_a5.scatter(Ns, [np.var(As) / np.mean(As) ** 2 for As in Ass], ec=colors.palette[0], fc="w", s=15, zorder=2)
    mean_vars = []
    var_means = []
    cv2s = []
    for N in Ns:
        mean_var = np.power(12 * rcls * rcls, -1) * (N + 1) * (N + 1) * (N + 2)
        var_mean = np.power(60 * rcls * rcls, -1) * (N + 1) * (N + 2) * (3 * N * N + 6 * N + 1) - np.power(
            36 * rcls * rcls, -1) * ((N + 1) * (N + 2)) ** 2
        mean = np.power(6 * rcls, -1) * (N + 1) * (N + 2)
        cv2 = (mean_var + var_mean) / (mean ** 2)
        mean_vars.append(mean_var / (mean ** 2))
        var_means.append(var_mean / (mean ** 2))
        cv2s.append(cv2)
    # ax_a5.plot(num_chas, [np.sqrt(MeanVarA + VarMeanA) for MeanVarA, VarMeanA in zip(Mean_of_VarA_over_n_theo, Var_of_MeanA_over_n_theo)], c=colors.palette[5], zorder=1)
    ax_a5.plot(Ns, cv2s, c=colors.palette[5], zorder=1)
    ax_a5.plot(Ns, mean_vars, c=colors.palette[5], ls=":")
    ax_a5.plot(Ns, var_means, c=colors.palette[5], ls="--")
    ax_a5.set_xticks(Ns)
    ax_a5.set_yticks([0.0, 0.5, 1.0])
    ax_a5.set_ylim([0, 1.5])

    r_ref = 5
    r_opn_single = 0.2
    # Subplot 1: x(t) over t
    ax_b1.set_xlabel("$t$ / s")
    ax_b1.set_ylabel("$x(t)$")
    data_no_ref = [(t, x, i) if x >= 0 else (t, 0, i) for (t, x, i) in data]
    data_plot = []
    for set1, set2 in zip(data_no_ref[:-1], data_no_ref[1:]):
        data_plot.append(set1)
        data_plot.append([set2[0], set1[1], set1[2]])

    ipis = get_ipis(data)
    starts_ipi = []
    stops_ipi = []
    is_ipi = False
    for elem in data:
        if elem[1] < 0 and not is_ipi:
            is_ipi = True
            t_ipi_start = elem[0]
            starts_ipi.append(t_ipi_start)
        if elem[1] > 0 and is_ipi:
            is_ipi = False
            t_ipi_stop = elem[0]
            stops_ipi.append(t_ipi_stop)

    for j in range(-2, 2):
        t_left = starts_ipi[puff_index1+ j] - t_start
        t_right = stops_ipi[puff_index1 +j] - t_start
        dt = t_right - t_left
        if j == -2:
            height = 0.5
            ax_b1.text((t_left + t_right) / 2, height + 1, "$...$", va="center", ha="center")
        elif j == -1:
            height = 0.5
            ax_b1.text((t_left + t_right) / 2, height +1, "$I_{i-1}$", va="center", ha="center")
        elif j == 0:
            height = 2.5
            ax_b1.text((t_left + t_right) / 2, height + 3, "$I_{i}$", va="center", ha="center")
        else:
            height = 0.5
            ax_b1.text((t_left + t_right) / 2, height +1, "$...$", va="center", ha="center")
        ax_b1.arrow(t_left, height, dt, 0, fc="k", length_includes_head=True, head_width=0.25, head_length=0.5,
                    lw=1.0, clip_on=False)
        ax_b1.arrow(t_right, height, -dt, 0, fc="k", length_includes_head=True, head_width=0.25,
                    head_length=0.5, lw=1.0, clip_on=False)

    ax_b1.plot(ts[t_start_index:t_stop_index] - t_start, xs[t_start_index:t_stop_index], lw=1, color=colors.palette[0])
    ax_b1.set_ylim([0, 7])
    ax_b1.set_xlim([0, t_stop - t_start])

    # Subplot 2: IPI density
    data_N1 = np.loadtxt(home + folder + "puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_K1_2.dat")
    ipis_N1 = get_ipis(data_N1)
    N = 2
    M = 3
    ropn = r_opn_single * N
    rref = r_ref
    mean = (M - 1) / rref + 1 / ropn
    var = (M - 1) * (1 / rref) ** 2 + (1 / ropn) ** 2
    dr = rref - ropn
    ts = np.linspace(0, max(ipis), 100)
    p_ipi = []
    for t in ts:
        c3 = np.power(rref / dr, 2) * (1 - np.exp(-dr * t) * (dr * t + 1))
        p3 = c3 * np.exp(-ropn * t)
        c4 = (1 / 2) * np.power(rref, 3) * (np.exp(-dr * t) * (-dr * t * (dr * t + 2) - 2) + 2) / np.power(dr, 3)
        p4 = c4 * np.exp(-ropn * t)
        p_ipi.append(ropn * p3)
    ax_b2.set_xlabel("$t$ / s")
    ax_b2.set_ylabel("$p_I(t; N=2)$")
    ax_b2.hist(ipis_N1, bins=50, color=colors.palette[0], density=True, label="$N=2$")
    ax_b2.plot(ts, p_ipi, color=colors.palette[5])
    ax_b2.set_ylim([0, 0.6])
    ax_b2.set_xlim([0, 10])

    # Subplot 2: IPI density
    data_N5 = np.loadtxt(home + folder + "puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_K1_5.dat")
    ipis_N5 = get_ipis(data_N5)
    N = 5
    M = 3
    ropn = r_opn_single * N
    rref = r_ref
    mean = (M - 1) / rref + 1 / ropn
    var = (M - 1) * (1 / rref) ** 2 + (1 / ropn) ** 2
    dr = rref - ropn
    ts = np.linspace(0, max(ipis), 100)
    p_ipi = []
    for t in ts:
        c3 = np.power(rref / dr, 2) * (1 - np.exp(-dr * t) * (dr * t + 1))
        p3 = c3 * np.exp(-ropn * t)
        c4 = (1 / 2) * np.power(rref, 3) * (np.exp(-dr * t) * (-dr * t * (dr * t + 2) - 2) + 2) / np.power(dr, 3)
        p4 = c4 * np.exp(-ropn * t)
        p_ipi.append(ropn * p3)
    ax_b3.set_xlabel("$t$ / s")
    ax_b3.set_ylabel("$p_I(t; N =5)$")
    ax_b3.hist(ipis_N5, bins=50, color=colors.palette[0], density=True, label="$N=5$")
    ax_b3.plot(ts, p_ipi, color=colors.palette[5])
    ax_b3.set_ylim([0, 0.6])
    ax_b3.set_xlim([0, 10])

    ts = np.linspace(0, 1, 101)
    Iss = []
    std_Is = []
    mean_Is = []
    CV_Is = []
    for N in Ns:
        file = f"puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_K1_{N:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        Is = get_ipis(data_n)
        Iss.append(Is)
        ropn = N * r_opn_single
        mean = (M - 1) / rref + 1 / ropn
        std = np.sqrt((M - 1) * (1 / rref) ** 2 + (1 / ropn) ** 2)
        mean_Is.append(mean)
        std_Is.append(std)
        CV_Is.append(std / mean)

    ax_b4.set_xlabel("$N$")
    ax_b4.set_xticks(Ns)
    ax_b4.set_ylabel(r"$\langle I \rangle^{-1}$ / s$^{-1}$")
    ax_b4.scatter(Ns, [1 / np.mean(Is) for Is in Iss], ec=colors.palette[0], fc="w", s=15, zorder=2)
    ax_b4.plot(Ns, [1 / mean_I for mean_I in mean_Is], c=colors.palette[5], zorder=1)
    ax_b4.set_yticks([0, 0.5, 1])
    ax_b4.set_ylim([0, 1.1])

    ax_b5.set_xlabel("$N$")
    ax_b5.set_xticks(Ns)
    ax_b5.set_ylabel(r"${CV}_I^2$")
    ax_b5.scatter(Ns, [np.var(Is) / np.mean(Is) ** 2 for Is in Iss], ec=colors.palette[0], fc="w", s=15, zorder=2)
    ax_b5.plot(Ns, [std_I ** 2 / mean_I ** 2 for std_I, mean_I in zip(std_Is, mean_Is)], c=colors.palette[5],
               zorder=1)
    ax_b5.set_yticks([0, 0.5, 1])
    ax_b5.set_ylim([0, 1.1])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig2b.pdf", transparent=True)
    plt.show()
