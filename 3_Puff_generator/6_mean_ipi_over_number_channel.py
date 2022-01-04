import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib import rcParams


import styles as st

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
    st.set_default_plot_style()
    home = os.path.expanduser("~")
    N = 10
    Ns = range(1, N)
    means = []
    for i in Ns:
        print(i)
        folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/ca_fix/"
        data = np.loadtxt(folder + "puff_markov_cafix0.33_ip1.00_tau1.00e+00_j1.00e+00_N1_{:d}.dat".format(i))
        data_tmp = []
        for set in data:
            if set[2] == 0:
                data_tmp.append(set)
        data = data_tmp
        ipis = ipi_distribution(data)
        means.append(np.mean(ipis))

    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax])
    ax.scatter([n+1 for n in Ns], [1/x for x in means])
    ax.set_xlabel("Cluster size")
    ax.set_ylabel("Mean event frequency")
    #plt.savefig(home + "/Data/Calcium/Plots/frequency_over_nr_cha.pdf", transparent=True)
    plt.show()