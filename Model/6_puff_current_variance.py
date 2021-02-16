import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

def autocorrelation(xs, dt):
    mean = np.mean(xs)
    xs = [x-mean for x in xs]
    correlations = []
    tmax = 1/2
    kmax = int(tmax/dt)
    ks = np.arange(0, kmax)
    for k in ks:
        print(k)
        corr = 0
        if k == 0:
            for x, y in zip(xs, xs):
                corr += x * y
        else:
            for x, y in zip(xs[:-k], xs[k:]):
                corr += x * y
        correlations.append(corr / len(xs[k:]))
    return [[dt * k for k in ks], correlations]


home = os.path.expanduser("~")
taus = np.logspace(-1, 2, 50)
jcas = np.logspace(-3, 0, 50)
CaFix = 0.95
Ncl = 10
NCha = 4
jca = 1
tau = 1
for i in range(1):
    data = np.loadtxt(home + "/CLionProjects/calcium-spikes-from-puff-phd/out/fixed calcium/spikes_fixedCa{:.2f}_adap_taua5.00e+01_e1.00e-01_tau{:.2e}_j{:.2e}_Ncha4_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(CaFix, tau, jca))
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

    dt0 = 0.1
    dts = []
    vars = []
    vars_dt = []
    factors = np.arange(1, 100)
    for f in factors:
        cg_jpuffs = coarse_grained_list(jpuffs_clear, f)
        vars.append(np.var(cg_jpuffs))
        vars_dt.append(np.var(cg_jpuffs)*dt0*f)
        dts.append(f*dt0)

    ax0.plot(dts, vars, label="$\sigma_{\Delta t}^2$")
    ax0.plot(dts, vars_dt, label= "$\sigma = \sigma_{\Delta t}^2 \Delta t$")

    # Cut data to only contain a certain puff site

    data = np.loadtxt(
        home + "/CLionProjects/calcium-spikes-from-puff-phd/out/fixed calcium/puff_fixedCa{:.2f}_adap_taua5.00e+01_e1.00e-01_tau1.00e+00_j1.00e+00_Ncha{:d}_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(
            CaFix, NCha))
    data_tmp = []
    for set in data:
        if set[2] == 2:
            data_tmp.append(set)
    data = data_tmp

    # Set all refractory states to 0
    for set in data[1:]:
        set[0] = set[0] - data[0][0]
    data[0][0] = 0
    for set in data:
        if set[1] < 0:
            set[1] = 0

    # Fill timesries of the cluster markov model
    ts = np.linspace(0, 5_000, 500_000)
    datas = []
    n = 1
    for t in ts:
        if t >= data[n][0]:
            n += 1
            if n > len(data) - 1:
                break
        datas.append(data[n - 1][1])

    dt = 0.01
    taus, corr = autocorrelation(datas, dt=dt)
    integral = 0
    deltats = np.logspace(-1, 1, 100)
    integrals = []
    for deltat in deltats:
        integral = 0
        for tau, c in zip(taus, corr):
            if tau < deltat:
                integral += c * (1 - tau / deltat) * dt
        integrals.append(integral)


    ax0.plot(deltats, [2*Ncl*x for x in integrals], label=r"$2 N \int_0^{\Delta t}d\tau C_{xx}(\tau)(1-\tau/\Delta t)$")
    ax0.set_xscale("log")

    ax0.legend()
    #ax0.set_ylim([0.01, 10])
    ax0.set_ylabel("$\sigma_{\Delta t}^2 , \sigma^2$")
    #ax0.set_yscale("log")
    ax0.set_xlabel("$\Delta t$")
    ax0.set_xscale("log")
    plt.savefig(home + "/Data/Calcium/Plots/variance_fixedCa{:.2f}_tau{:.2e}_j{:.2e}_nClu10_nCha4_adap.pdf".format(CaFix, tau, jca), transparent=True)
    plt.show()
