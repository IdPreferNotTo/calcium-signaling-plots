import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

def ipi_distribution(data):
    minimum = min([x[1] for x in data])
    isis = []
    t_tmp = 0
    for set in data: #set  = [time, state]
        if set[1] == minimum:
            isis.append(set[0] - t_tmp)
            t_tmp = set[0]
        #elif(set[1] > 0):
        #    isis.append(set[0] - t_tmp)
    return isis


if __name__ == "__main__":
    set_default_plot_style()
    home = os.path.expanduser("~")
    N = 10
    Ns = range(N)
    means = []
    for i in Ns:
        print(i)
        data = np.loadtxt(
            home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N10_{:d}.dat".format(i))
        data_tmp = []
        for set in data:
            if set[2] == 3:
                data_tmp.append(set)
        data = data_tmp
        ipis = ipi_distribution(data)
        means.append(np.mean(ipis))

    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    remove_top_right_axis([ax])
    ax.scatter([n+1 for n in Ns], [1/x for x in means])
    ax.set_xlabel("Cluster size")
    ax.set_ylabel("Mean event frequency")
    plt.savefig(home + "/Data/Calcium/Plots/frequency_over_nr_cha.pdf", transparent=True)
    plt.show()