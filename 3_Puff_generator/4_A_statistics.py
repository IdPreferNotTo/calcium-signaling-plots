import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
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

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4.5, 3))
    grids = gridspec.GridSpec(2, 2)
    ax0 = fig.add_subplot(grids[0, 0])
    ax1 = fig.add_subplot(grids[0, 1])
    ax2 = fig.add_subplot(grids[1, 0])
    ax3 = fig.add_subplot(grids[1, 1])

    axis = [ax0, ax1, ax2, ax3]
    st.remove_top_right_axis(axis)
    home = os.path.expanduser("~")

    # Subplot 1: x(t) over t
    cafix = 0.2
    ax0.set_xlabel("$t$ / s")
    ax0.set_ylabel("$x(t)$")
    folder = "/Data/calcium_spikes_markov/ca_fix/"
    data = np.loadtxt(home + folder + f"puff_markov_cafix{cafix:.2f}_ip1.00_tau1.00e+00_j1.00e+00_K1_5.dat")
    data_tmp = []

    data_no_ref = [(t, x, i) if x>= 0 else (t, 0, i) for (t, x, i, sx) in data]
    data2 = []
    for set1, set2 in zip(data_no_ref[:-1], data_no_ref[1:]):
        data2.append(set1)
        data2.append([set2[0], set1[1], set1[2]])
    ts, xs, idxs = np.transpose(data2)

    starts_ipi =  []
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

    ax0.plot(ts[:500] - ts[0], xs[:500], lw=1, color=st.colors[0])
    ax0.set_xlim([t_left-0.15, t_right+0.15])
    ax0.fill_between(ts[:500] - ts[0], xs[:500], color=st.colors[0], alpha=0.5)

    ax1.set_ylim([0, 1/5 + 0.15])
    ax1.set_xlabel(r"$n_0$")
    ax1.set_ylabel(r"$p(n_0)$")
    ns_open = get_ns_open(data)
    ns_over_i = [0]*5
    for n in ns_open:
        ns_over_i[n-1] += 1
    for i, n in enumerate(ns_over_i):
        ax1.bar(i+1, n / len(ns_open), 0.8, align="center", color=st.colors[0])
    ax1.set_xticks([1, 2, 3, 4, 5])
    ax1.set_yticks([0, 0.1, 0.2, 0.3])
    #ax1.plot([0, 1, 1, 2, 3, 4, 5, 5, 6], [0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0], lw=1, color=st.colors[2])

    ax2.set_xlabel("$A$")
    ax2.set_ylabel("$p(A)$")
    As = get_As(data_no_ref)
    ax2.hist(As, bins=25, color=st.colors[0], density=True)


    Ass = []
    num_chas = np.arange(1, 7)
    for num_cha in num_chas:
        file = f"puff_markov_cafix{cafix:.2f}_ip1.00_tau1.00e+00_j1.00e+00_K1_{num_cha:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        As = get_As(data_n)
        Ass.append(As)

    ax3.set_xlabel("$N$")
    ax3.set_ylabel(r"$\langle A \rangle$")
    ax3.scatter(num_chas, [np.mean(As) for As in Ass], ec=st.colors[0], fc="w", s=15, zorder=2)
    Ass_theo = [(1/50)*(n +1)*(n+2)/6 for n in num_chas]
    ax3.plot(num_chas, Ass_theo, lw=1, c=st.colors[2], zorder=1)
    ax3.set_xticks([1, 2, 3, 4, 5, 6])
    plt.show()
