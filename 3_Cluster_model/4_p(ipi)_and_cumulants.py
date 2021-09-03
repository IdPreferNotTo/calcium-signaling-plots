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


def get_ipis_from_single_cluster(data, nr_cluster):
    # data = [t, state, idx]

    # Cut data to contain data from a a single cluster only
    single_cluster = []
    for elem in data:
        if elem[2] == nr_cluster:
            single_cluster.append(elem)

    # Turn data into interpuff intervals
    t0 = single_cluster[0][0]
    ipis= []
    for elem in single_cluster:
        if elem[1] == 0:
            ipi = elem[0] - t0
            ipis.append(ipi)
            t0 = elem[0]
    return ipis


if __name__ == "__main__":
    set_default_plot_style()
    home = os.path.expanduser("~")
    fig = plt.figure(tight_layout=True, figsize=(4, 6*(3/4)))

    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0])
    #ax1in = inset_axes(ax1, width=1.5, height=0.8)
    ax2 = fig.add_subplot(gs[1])
    remove_top_right_axis([ax1, ax2])

    num_cls = 4
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/"
    file = "puff_markov_cafix0.33_ip1.00_tau1.00e-01_j1.00e-03_N1_Ncls{:d}.dat".format(num_cls)
    data = np.loadtxt(home + folder + file)

    ipis = get_ipis_from_single_cluster(data, 0)
    mean_ipi = np.mean(ipis)
    CV = np.std(ipis)/np.mean(ipis)
    ax1.hist(ipis, bins=25, density=True, alpha=0.7, label=rf"$n_{{\rm cls}} = {num_cls}$" + "\n" + f"$\mu(I) = {mean_ipi:.2f}$" + "\n" + f"$C_V(I) = {CV:.2f}$")
    m=4
    n=4

    mean_ibi = 10
    ratio = 0.1
    ropn = (1. + ratio * (num_cls - 1.)) / mean_ibi
    rref = (1. + ratio * (num_cls - 1.)) / (ratio * mean_ibi)
    mean = (m-1)/rref + 1/ropn
    var = (m-1)*(1/rref)**2 + (1/ropn)**2
    dr = rref - ropn

    ts = np.linspace(0, max(ipis), 100)
    p_ipi = []
    if num_cls == 1:
        for t in ts:
            p1 = np.exp(-ropn*t)
            p_ipi.append(ropn*p1)
    elif num_cls == 2:
        for t in ts:
            p2 = (rref/dr)*(1 - np.exp(-dr*t))*np.exp(-ropn*t)
            p_ipi.append(ropn*p2)
    elif num_cls == 3:
        for t in ts:
            p3  = np.power((rref/dr), 2)*(1 - np.exp(-dr*t)*(dr*t + 1))*np.exp(-ropn*t)
            p_ipi.append(ropn*p3)
    elif num_cls == 4:
        for t in ts:
            c4 = (1/2)*np.power(rref,3)*(np.exp(-dr*t)*(-dr*t*(dr*t +2) -2) +2)/np.power(dr, 3)
            p4 = c4*np.exp(-ropn*t)
            p_ipi.append(ropn*p4)

    ax1.plot(ts, p_ipi, c="k")
    print(mean, var)
    print(np.mean(ipis), np.var(ipis))
    ax1.legend()
    ax1.set_xlabel("$I$")
    ax1.set_ylabel("$p(I)$")

    for num_cha in range(2, 7):
        folder = "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/"
        file = "puff_markov_cafix0.33_ip1.00_tau1.00e-01_j1.00e-03_N1_{:d}.dat".format(num_cha)
        data = np.loadtxt(home + folder + file)
        ipis = get_ipis_from_single_cluster(data, 0)
        r = 0.1
        CV2 = (1. + (m-1) * r**2)/(1. + (m-1)*r)**2
        CV = np.sqrt(CV2)
        xs = np.linspace(0, 5, 100)
        ax2.plot(xs, [CV*x for x in xs], c="C7", zorder=4)
        ax2.plot(xs, xs, ls="--", c="C7", zorder=4)
        for i in range(5):
            if 100*(i+1) > len(ipis):
                break
            else:
                ipis_tmp = ipis[100*i:100*(i+1)]
                ax2.scatter(np.mean(ipis_tmp), np.std(ipis_tmp), c="C0")
    ax2.set_xlim([0, 5])
    ax2.set_ylim([0, 5])
    ax2.set_xlabel("$\mu$")
    ax2.set_ylabel("$\sigma$")
    plt.savefig(
        home + "/Data/Calcium/Plots/ipi_cv.pdf", transparent=True)
    plt.show()