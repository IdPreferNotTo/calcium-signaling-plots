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

    def coarse_grained_list(list, factor):
        coarse_grained_list = []
        max = int(len(list)/factor)
        for i in range(max-1):
            coarse = np.mean(list[factor*i:factor*(i+1)])
            coarse_grained_list.append(coarse)
        return coarse_grained_list

    jpuffs_dt01 = jpuffs_clear
    jpuffs_dt1 = coarse_grained_list(jpuffs_clear, 10)
    jpuffs_dt10 = coarse_grained_list(jpuffs_clear, 100)

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(2, 2)

    ax_ts1 = fig.add_subplot(gs[0, 0])
    ax_hist1 = fig.add_subplot(gs[0, 1])
    ax_ts2 = fig.add_subplot(gs[1, 0])
    ax_hist2 = fig.add_subplot(gs[1, 1])
    plt_start = 0
    plt_stop = plt_start + 500

    mean = np.mean(jpuffs_clear)
    ax_ts1.set_ylabel(r"J$_{\rm Puff}$")
    ax_ts1.plot(ts_clear[plt_start:plt_stop], jpuffs_clear[plt_start:plt_stop], c="C2",label="$\mu_0 = {:.2f}$".format(mean_jpuffs_clear[1]))
    ax_ts1.plot(ts_clear[plt_start:plt_stop], mean_jpuffs_clear[plt_start:plt_stop], c="k", label="$\mu = {:.2f}$".format(mean))
    ax_ts1.set_xlabel("$t$")
    ax_ts1.set_xlim([0, 50])
    ax_ts1.legend()

    mean = np.mean(jpuffs_clear)
    ax_ts2.set_ylabel(r"J$_{\rm Puff}$")
    ts_clear_dt1 = [10*t for t in ts_clear]
    ax_ts2.plot(ts_clear_dt1[plt_start:plt_stop], jpuffs_dt1[plt_start:plt_stop], c="C3",label="$\mu_0 = {:.2f}$".format(mean_jpuffs_clear[1]))
    ax_ts2.plot(ts_clear_dt1[plt_start:plt_stop], mean_jpuffs_clear[plt_start:plt_stop], c="k", label="$\mu = {:.2f}$".format(mean))
    ax_ts2.set_xlabel("$t$")
    ax_ts2.set_xlim([0, 50])
    ax_ts2.legend()

    std01 = np.std(jpuffs_dt01)
    std1 = np.std(jpuffs_dt1)
    std10 = np.std(jpuffs_dt10)
    var01 = np.var(jpuffs_dt01)
    var1 = np.var(jpuffs_dt1)
    var10 = np.var(jpuffs_dt10)

    gauss001 = np.linspace(mean - 3*std01, mean + 3*std01, 100)
    gauss01 = np.linspace(mean - 3*std1, mean + 3*std1, 100)
    def gauss_dist(xs, mean, std):
        gauss_dist = []
        for x in xs:
            gauss = 1/np.sqrt(2*np.pi*(std**2)) * np.exp(-((x - mean)**2)/(2*std**2))
            gauss_dist.append(gauss)
        return gauss_dist

    ax_hist1.plot(gauss001, gauss_dist(gauss001, mean, std01), c="C2")
    ax_hist1.hist(jpuffs_dt01, bins = 20, alpha = 0.7, color="C2", density=True, label="$\Delta t$ = {:.2f} \n $\sigma_1$ = {:.3f}".format(0.1, std01))
    ax_hist1.set_xlabel("$y_{0.01}$")
    ax_hist1.set_ylabel("$P(y)$")
    ax_hist1.legend()



    ax_hist2.plot(gauss01, gauss_dist(gauss01, mean, std1), c="C3")
    ax_hist2.hist(jpuffs_dt1, bins=20, alpha=0.7, color="C3", density=True,
             label="$\Delta t$ = {:.2f} \n $\sigma_2$ = {:.3f}".format(1.0, std1))
    ax_hist2.set_xlabel("$y_1$")
    ax_hist2.set_ylabel("$P(y)$")
    ax_hist2.legend()




    plt.savefig(home + "/Data/Calcium/Plots/correlations_fixedCa{:.2f}_tau{:.2e}_j{:.2e}_nClu10_nCha4_adap.pdf".format(fixedCa, tau, jca), transparent=True)
    plt.show()
