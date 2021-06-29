import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def ipi_distribution(data):
    minimum = min([x[1] for x in data])
    isis = []
    t_tmp = 0
    for set in data: #set  = [time, state]
        if set[1] == minimum:
            t_tmp = set[0]
        elif(set[1] > 0):
            isis.append(set[0] - t_tmp)
    return isis


if __name__ == "__main__":
    home = os.path.expanduser("~")
    fig = plt.figure(tight_layout=True, figsize=(9 / 2, 6))
    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    data = np.loadtxt(
        home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/puff_markov_cafix0.33_tau1.00e+00_j1.00e+00_N10_0.dat")

    data_tmp = []
    for set in data:
        if set[2] == 3:
            data_tmp.append(set)
    data = data_tmp
    ipis = ipi_distribution(data)
    ax1.hist(ipis, bins=20, density=True,
                label="$n = 4$ \n $\mu(I) = {:.2f}$ \n $C_V(I) = {:.2f}$".format(np.mean(ipis), np.std(ipis) / np.mean(ipis)))
    m=4
    n=4
    rref = 1.30 * n
    ropn = 0.13 * n
    mean = (m-1)/rref + 1/ropn
    var = (m-1)*(1/rref)**2 + (1/ropn)**2
    print(mean, var)
    print(np.mean(ipis), np.var(ipis))
    ax1.legend()
    ax1.set_xlabel("$I$")
    ax1.set_ylabel("$p(I)$")


    means = []
    stds = []
    N = 7
    for i in range(1, N):
        data = np.loadtxt(
            home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/puff_markov_cafix0.33_tau1.00e+00_j1.00e+00_N10_0.dat")

        data_tmp = []
        for set in data:
            if set[2] == 3:
                data_tmp.append(set)
        data = data_tmp
        ipis = ipi_distribution(data)
        for _ in range(5):
            rng = np.random.randint(500)
            mean = np.mean(ipis[rng:rng+100])
            std = np.std(ipis[rng:rng+100])
            means.append(mean)
            stds.append(std)
    ax2.scatter(means, stds, zorder=1)
    r = 0.1
    CV2 = (1. + (m-1) * r**2)/(1. + (m-1)*r)**2
    CV = np.sqrt(CV2)
    xs = np.linspace(0, 5, 100)
    ax2.plot(xs, [CV*x for x in xs], c="C7", zorder=4, label="slope = {:.2f}".format(CV))
    ax2.plot(xs, xs, ls="--", c="C7", zorder=4)
    ax2.set_xlim([0, 5])
    ax2.set_ylim([0, 5])
    ax2.set_xlabel("$\mu$")
    ax2.set_ylabel("$\sigma$")
    ax2.legend()
    plt.savefig(
        home + "/Data/Calcium/Plots/ipi_cv.pdf", transparent=True)
    plt.show()