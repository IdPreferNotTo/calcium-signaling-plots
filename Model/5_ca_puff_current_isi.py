import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


home = os.path.expanduser("~")
taus = np.logspace(-1, 2, 50)
jcas = np.logspace(-3, 0, 50)
Ncl = 10
Nch = 4
for i in range(3):
    tau = [taus[30], taus[35], taus[40]][i] #taus[45]
    jca = jcas[25] #[jcas[15], jcas[25], jcas[35]][i]
    data = np.loadtxt(home + "/CLionProjects/calcium-spikes-from-puff-phd/out/Data/spikes_adap_taua5.00e+01_e1.00e-01_tau{:.2e}_j{:.2e}_Ncha4_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(tau, jca))
    ts, cas, jpuffs, adaps = np.transpose(data)

    t_spike = 0
    ISI = []
    ts_abs = []
    t_tmp = 0
    for t, Ca in zip(ts, cas):
        if Ca == 1:
            ISI.append(t - t_tmp)
            t_tmp = t
    mean_ISI = np.mean(ISI)
    cv_ISI = np.var(ISI) / (mean_ISI ** 2)

    ts_clear = []
    cas_clear = []
    jpuffs_clear = []
    mean_jpuffs_clear = []
    for t, ca, jpuff, adap in zip(ts, cas, jpuffs, adaps):
        if ca!=1:
            ts_clear.append(t)
            cas_clear.append(ca)
            jpuffs_clear.append(jpuff)
            mean_jpuff = jca*(Ncl*(Nch + 2)/3)*np.power(1 + 2*np.power(Nch*(Nch+1), -1, dtype=float)*np.power(0.33/ca, 3)*(1+ca**3)/(1+0.33**3)*10/0.02, -1)
            mean_jpuffs_clear.append(mean_jpuff)


    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(7, 1)

    ax0 = fig.add_subplot(gs[0:2])
    ax1 = fig.add_subplot(gs[2:5])
    ax2 = fig.add_subplot(gs[5:7])
    plt_start = 0
    plt_stop = plt_start + int(10*mean_ISI*10)

    ax0.set_ylabel(r"J$_{\rm Puff}$")
    ax0.plot(ts_clear[plt_start:plt_stop], jpuffs_clear[plt_start:plt_stop], c="C2")
    ax0.plot(ts_clear[plt_start:plt_stop], mean_jpuffs_clear[plt_start:plt_stop], c="k")

    ax1.set_ylabel("[Ca$^{2+}$]")
    ax1.set_xlabel(r"$t$")

    ax1.plot(ts[plt_start:plt_stop], cas[plt_start:plt_stop],
             label="$T^* = {:.2f}$".format(mean_ISI) + "\n $C_V = {:.2f}$".format(np.sqrt(cv_ISI)))
    ax1.legend()

    ax2.set_ylabel(r"$\rm{ISI} / \langle \rm{ISI} \rangle$")
    ax2.scatter(range(20), ISI[:20]/mean_ISI)

    fig.align_ylabels([ax0, ax1, ax2])

    plt.savefig(home + "/Data/Calcium/Plots/Ca_tau{:.2e}_j{:.2e}_nClu10_nCha4_adap.pdf".format(tau, jca), transparent=True)
    plt.show()
