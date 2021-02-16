import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


home = os.path.expanduser("~")
taus = np.logspace(-1, 2, 50)
jcas = np.logspace(-3, 0, 50)
fixedCa = 0.35
Ncl = 10
Nch = 4
jca = 1
tau = 1
for i in range(1):
    data = np.loadtxt(home + "/CLionProjects/calcium-spikes-from-puff-phd/out/fixed calcium/spikes_fixedCa{:.2f}_adap_taua5.00e+01_e1.00e-01_tau{:.2e}_j{:.2e}_Ncha4_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(fixedCa, tau, jca))
    ts, cas, jleaks, jpuffs = np.transpose(data)

    ts_clear = []
    cas_clear = []
    jpuffs_clear = []
    mean_jpuffs_clear = []
    for t, ca, jpuff in zip(ts, cas, jpuffs):
        if ca!=1:
            ts_clear.append(t)
            cas_clear.append(ca)
            jpuffs_clear.append(jpuff)
            Popen = np.power(1 + 2*np.power(Nch*(Nch+1), -1, dtype=float)*np.power(0.33/ca, 3)*(1+ca**3)/(1+0.33**3)*10/0.02, -1)
            mean_jpuff = jca*(Ncl*(Nch + 2)/3)*Popen
            mean_jpuffs_clear.append(mean_jpuff)

    std_jpuff = Popen*(Nch+1)*(Nch+2)/3
    print(std_jpuff/np.sqrt(0.1))

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(5, 2)

    ax2 = fig.add_subplot(gs[0:2, 1])
    ax1 = fig.add_subplot(gs[0:2, 0])
    ax0 = fig.add_subplot(gs[2:5, :2])
    plt_start = 0
    plt_stop = plt_start + 1000

    ax0.set_ylabel(r"J$_{\rm Puff}$")
    mean = np.mean(jpuffs_clear)
    ax0.plot(ts_clear[plt_start:plt_stop], jpuffs_clear[plt_start:plt_stop], c="C2",label="$\mu_0 = {:.2f}$".format(mean_jpuffs_clear[1]))
    ax0.plot(ts_clear[plt_start:plt_stop], mean_jpuffs_clear[plt_start:plt_stop], c="k", label="$\mu = {:.2f}$".format(mean))
    ax0.legend()

    ax1.set_ylabel("[Ca$^{2+}$]")
    ax1.set_xlabel(r"$t$")

    ax1.plot(ts[plt_start:plt_stop], cas[plt_start:plt_stop])
    ax1.legend()

    def coarse_grained_list(list, factor):
        coarse_grained_list = []
        max = int(len(list)/factor)
        for i in range(max-1):
            coarse = np.mean(list[factor*i:factor*(i+1)])
            coarse_grained_list.append(coarse)
        return coarse_grained_list

    jpuffs_dt001 = jpuffs_clear
    jpuffs_dt01 = coarse_grained_list(jpuffs_clear, 10)
    jpuffs_dt1 = coarse_grained_list(jpuffs_clear, 100)
    std001 = np.std(jpuffs_dt001)
    std01 = np.std(jpuffs_dt01)
    gauss001 = np.linspace(mean - 3*std001, mean + 3*std001, 100)
    gauss01 = np.linspace(mean - 3*std01, mean + 3*std01, 100)
    def gauss_dist(xs, mean, std):
        gauss_dist = []
        for x in xs:
            gauss = 1/np.sqrt(2*np.pi*(std**2)) * np.exp(-((x - mean)**2)/(2*std**2))
            gauss_dist.append(gauss)
        return gauss_dist

    ax2.plot(gauss001, gauss_dist(gauss001, mean, std001), c="C2")
    ax2.plot(gauss01, gauss_dist(gauss01, mean, std01), c="C3")
    ax2.hist(jpuffs_dt001, bins = 20, alpha = 0.7, color="C2", density=True, label="dt = {:.2f} \n $\sigma$ = {:.2f}".format(0.01, std001))
    ax2.hist(jpuffs_dt01, bins=20, alpha=0.7, color="C3", density=True, label="dt = {:.2f} \n $\sigma$ = {:.2f}".format(0.1, std01))
    #ax2.hist(jpuffs_dt1, bins=20, density=True, label="dt = {:.2f}".format(1.0))
    ax2.legend()

    fig.align_ylabels([ax0, ax1])

    plt.savefig(home + "/Data/Calcium/Plots/fixedCa{:.2f}_tau{:.2e}_j{:.2e}_nClu10_nCha4_adap.pdf".format(fixedCa, tau, jca), transparent=True)
    plt.show()
