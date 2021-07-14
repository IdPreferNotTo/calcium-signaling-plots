import os
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
from matplotlib import rcParams


def set_default_plot_style():
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = 'Computer Modern Roman'
        rc('text', usetex=True)


def remove_top_right_axis(axis):
        for ax in axis:
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)


if __name__ == "__main__":
    set_default_plot_style()
    home = os.path.expanduser("~")
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    axin = inset_axes(ax, width=1.5, height=0.8)
    remove_top_right_axis([ax])
    n = 4
    m = 4
    data = np.loadtxt(home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/puff_markov_cafix0.33_tau1.00e+00_j1.00e+00_N10_0.dat")

    data_tmp = []
    for set in data:
        if set[2] == 3:
            data_tmp.append(set)
    data = data_tmp
    data2 = []
    for set1, set2 in zip(data[:-1], data[1:]):
        data2.append(set1)
        data2.append([set2[0], set1[1], set1[2]])

    ts, xs, idxs = np.transpose(data2)
    axin.plot(ts, xs)
    axin.set_xlim([55, 60])
    axin.tick_params(direction="in")
    axin.set_xticklabels([])
    axin.set_yticklabels([])

    ys =[]
    for x in xs:
        if x > 0:
            ys.append(x)
        else:
            ys.append(0)
    ax.plot(ts, ys)
    rect = patches.Rectangle((55, -0.1), 5, 4.2, linewidth=1, edgecolor='k', facecolor='none')
    ax.add_patch(rect)

    ax.set_ylim([-0.2, 5])
    ax.set_xlim([50, 100])
    ax.set_xlabel("time $t$")
    ax.set_ylabel("$x(t)$")
    plt.savefig(home + "/Data/Calcium/Plots/cluster_timeseries.pdf", transparent=True)
    plt.show()