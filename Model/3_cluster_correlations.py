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
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(2, 2)
ax = fig.add_subplot(gs[0, :])
ax_corr = fig.add_subplot(gs[1, 0])
ax_comp = fig.add_subplot(gs[1, 1])
ax.set_xlabel("$t$")
ax.set_ylabel("$x(t)$")

#double kCaRest = 0.33; // Ca is given in terms of Kca (r_open = r_max (Ca^n / Ca^n + Kca^n) (Ip3^m / Ip3^m + Kip3^m))
CaRest = 0.33
CaFix = 0.35
Ip3 = 1.
NCha = 4
NClosedStates = 4
m = NClosedStates - 2
meanI = 10 # for a single channel

data = np.loadtxt(home + "/CLionProjects/calcium-spikes-from-puff-phd/out/fixed calcium/puff_fixedCa{:.2f}_adap_taua5.00e+01_e1.00e-01_tau1.00e+00_j1.00e+00_Ncha{:d}_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(CaFix, NCha))
data_tmp = []
for set in data:
    if set[2] == 5:
        data_tmp.append(set)
data = data_tmp
data2 = []
for set1, set2 in zip(data[:-1], data[1:]):
    data2.append(set1)
    data2.append([set2[0], set1[1], set1[2]])

x, t, idx = np.transpose(data2)
ax.plot(x, t)
ax.set_xlim([50, 100])
ax.set_yticks([-3, -2, -1, 0, 1, 2, 3, 4])
ax.set_yticklabels(["$0_4$", "$0_3$", "$0_2$", "$0_1$", "$1$", "$2$", "$3$", "$4$"])

ts = np.linspace(0, 5_000, 500_000)
datas = []

for set in data[1:]:
    set[0] = set[0] - data[0][0]
data[0][0] = 0
for set in data:
    if set[1] < 0:
        set[1] = 0
n = 1
for t in ts:
    if t >= data[n][0]:
        n += 1
        if n > len(data)-1:
            break
    datas.append(data[n-1][1])

dt=0.01
taus, corr = autocorrelation(datas, dt=dt)
integral = 0
deltats =np.logspace(-2, 0, 100)
integrals = []
for deltat in deltats:
    integral = 0
    for tau, c in zip(taus,corr):
        if tau < deltat:
            integral += c*(1 - tau/deltat)*dt
    integrals.append(integral)
ax_comp.plot(deltats, integrals)
ax_comp.axhline(y=sum(corr)*dt, c="k")
#ax_comp.set_ylim([0, 0.005])
ax_comp.set_xscale("log")

ax_corr.plot(taus, corr, label=r"$\int d\tau C_{{xx}}(\tau) = {:.2e}$".format(sum(corr)*dt))
ax_corr.set_xlabel(r"$\tau$")
ax_corr.set_ylabel(r"$C_{xx}(\tau)$")
ax_corr.fill_between(taus, corr, 0, facecolor="C0", alpha=0.5, zorder=2)
ax_corr.legend()


plt.savefig(home + "/Data/Calcium/Plots/single_cluster_correlations_fixedCa{:.2f}_meanI{:.1f}_ncls{:d}_ncha{:d}.pdf".format(CaFix, meanI, NClosedStates, NCha), transparent=True)
plt.show()