import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

def moments(xs, k):
    moment = 0
    mu = np.mean(xs)
    for x in xs:
        moment += (x - mu)**k
    return moment/len(xs)

home = os.path.expanduser("~")
fixedCa = 0.33
Ncl = 10
Nch = 4
jca = 1
tau = 1
for i in range(1):
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/"
    file = "ca_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N10_0.dat".format(fixedCa)
    data = np.loadtxt(home + folder + file)

    ts, cas, jpuffs, adaps = np.transpose(data)
    ts_clear = []
    cas_clear = []
    jpuffs_clear = []
    mean_jpuffs_clear = []
    for t, ca, jpuff in zip(ts, cas, jpuffs):
        if ca!=1:
            ts_clear.append(t)
            cas_clear.append(ca)
            jpuffs_clear.append(jpuff)
    mean = np.mean(jpuffs_clear)

    def coarse_grained_list(list, factor):
        coarse_grained_list = []
        max = int(len(list)/factor)
        for i in range(max-1):
            coarse = np.mean(list[factor*i:factor*(i+1)])
            coarse_grained_list.append(coarse)
        return coarse_grained_list

    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax0 = fig.add_subplot(gs[0, 0])

    dt0 = 0.001
    dts = []
    mu2s = []
    mu3s = []
    mu4s = []

    factors = np.arange(1, 100)
    for f in factors:
        cg_jpuffs = coarse_grained_list(jpuffs_clear, f)
        mu2 = moments(cg_jpuffs, 2)
        mu2s.append(mu2)
        mu3 = moments(cg_jpuffs, 3)
        mu3s.append(mu3)
        mu4 = moments(cg_jpuffs, 4)
        mu4s.append(mu4)
        dts.append(f*dt0)

    ax0.plot(dts, mu2s, label=r"$\gamma_2 = \langle (x - \mu)^2  \rangle$")
    ax0.plot(dts, mu3s, label=r"$\gamma_3 = \langle (x - \mu)^3  \rangle$")
    ax0.plot(dts, mu4s, label=r"$\gamma_4 = \langle (x - \mu)^4  \rangle$")
    ax0.legend()
    #ax0.set_ylim([0.001, 1])
    ax0.set_ylabel("$\mu_x$")
    ax0.set_yscale("log")
    ax0.set_xlabel("$\Delta t$")
    ax0.set_xscale("log")
    plt.savefig(home + "/Data/Calcium/Plots/markov_jpuff_moments_over_deltat_static.pdf".format(fixedCa, tau, jca), transparent=True)
    plt.show()
