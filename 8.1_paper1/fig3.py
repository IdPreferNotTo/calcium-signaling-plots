import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches
import os

import styles as st

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
    fig = plt.figure(tight_layout=True, figsize=(4.5, 3))
    grids = gridspec.GridSpec(2, 2)
    ax0 = fig.add_subplot(grids[0, :])
    ax1 = fig.add_subplot(grids[1, 0])
    ax2 =  fig.add_subplot(grids[1, 1])
    axis = [ax0, ax1, ax2]
    st.remove_top_right_axis(axis)
    home = os.path.expanduser("~")

    # Subplot 1: x(t) over t
    ax0.set_xlabel("$t$ / s")
    ax0.set_ylabel("$x(t)$")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/"
    data = np.loadtxt(home + folder + "puff_markov_cafix0.33_ip1.00_taua1.00e+00_ampa1.00e+00_tau1.00e+00_j1.00e+00_N1_5.dat")

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

    t_left = starts_ipi[8] - ts[0]
    t_right = stops_ipi[8] - ts[0]
    dt = t_right - t_left
    ax0.arrow(t_left, 2.8, dt, 0, fc="k", length_includes_head=True, head_width=0.5, head_length=1.0,
              lw=0.5, clip_on=False)
    ax0.arrow(t_right, 2.8, -dt, 0, fc="k", length_includes_head=True, head_width=0.5,
              head_length=1.0, lw=0.5, clip_on=False)
    ax0.text((t_left + t_right) / 2, 3.5, "$I_i$", va="center", ha="center")

    ax0.plot(ts[:500] - ts[0], xs[:500], lw=1, color=st.colors[0])
    ax0.set_xlim([0, 25])


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
    ax1.set_xlabel("$I$ / s")
    ax1.set_ylabel("$P(I)$")
    ax1.hist(ipis, bins=50, color=st.colors[0], density=True)
    ax1.plot(ts, p_ipi, color=st.colors[2])
    ax1.set_xlim([0, 10])

    ax2.set_xlabel(r"$\langle I \rangle$ / s")
    ax2.set_ylabel(r"$\sqrt{\langle \Delta I^2 \rangle}$ / s")
    Ass = []
    num_chas = np.arange(1, 7)
    for num_cha in num_chas:
        file = f"puff_markov_cafix0.33_ip1.00_taua1.00e+00_ampa1.00e+00_tau1.00e+00_j1.00e+00_N1_{num_cha:d}.dat"
        data_n = np.loadtxt(home + folder + file)
        ipis = get_ipis(data_n)
        r = 0.1
        CV2 = (1. + (num_cls-1) * r**2)/(1. + (num_cls-1)*r)**2
        CV = np.sqrt(CV2)
        xs = np.linspace(0, 6, 100)
        ax2.plot(xs, [CV*x for x in xs], lw=1, c=st.colors[2], zorder=1)
        ax2.plot(xs, xs, lw=1, ls=":", c="k", zorder=1)
        for i in range(5):
            if 100*(i+1) > len(ipis):
                break
            else:
                ipis_tmp = ipis[100*i:100*(i+1)]
                ax2.scatter(np.mean(ipis_tmp), np.std(ipis_tmp), ec=st.colors[0], fc="w", s=15, zorder=2)
    ax2.set_xlim([0, 6])
    ax2.set_ylim([0, 6])
    style = "Simple, head_width=3, head_length=5"
    kw = dict(arrowstyle=style, color="k")
    arrow = patches.FancyArrowPatch((4.2, 4.2), (5, 0),
                             connectionstyle="arc3,rad=-0.25", zorder=5, **kw)
    ax2.add_patch(arrow)
    ax2.text(5.2, 1.5, "$m$")
    ts = np.linspace(0, 1, 101)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig3.pdf",
                transparent=True)
    plt.show()
