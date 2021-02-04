import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

def ipi_distribution(data):
    minimum = min([x[1] for x in data])
    print(minimum)
    isis = []
    t_tmp = 0
    for set in data: #set  = [time, state]
        if set[1] == minimum:
            t_tmp = set[0]
        elif(set[1] > 0):
            isis.append(set[0] - t_tmp)
    return isis

def amp_distribution(data):
    amps = []
    for set1, set2 in zip(data[:-1], data[1:]):
        if set1[1] == 0:
            amps.append(set2[1])
    return amps

home = os.path.expanduser("~")
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(2, 2)
ax = fig.add_subplot(gs[0,:])
ax_isi = fig.add_subplot(gs[1,0])
ax_amp = fig.add_subplot(gs[1,1])
ax.set_xlabel("$t$")
ax.set_ylabel("$x(t)$")
ax_isi.set_xlabel("$I$")
ax_isi.set_ylabel("$p(I)$")
ax_amp.set_xlabel("Amp.")
ax_amp.set_ylabel(r"$p(Amp.)$")

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
    if set[2] == 3:
        data_tmp.append(set)
data = data_tmp
data2 = []
for set1, set2 in zip(data[:-1], data[1:]):
    data2.append(set1)
    data2.append([set2[0], set1[1], set1[2]])

x, t, idx = np.transpose(data2)
ax.plot(x, t)
ax.set_xlim([50, 100])

ts  = np.linspace(0, 20, 200)
ps = []
rrefsingle = 1.3*(CaFix/CaRest)**3 *(1 + CaRest**3)/(1 + CaFix**3)
ropensingle = 0.13*(CaFix/CaRest)**3 *(1 + CaRest**3)/(1 + CaFix**3)
rref = ropensingle+NCha
ropen = ropensingle*NCha
x = rref - ropen
for t in ts:
    p = ropen*0.5*(rref**(m+1))*(np.exp(-rref*t)*(-(t*x)**m - 2*t*x -2)  + 2*np.exp(-ropen*t) )/(x**(m+1))
    ps.append(p)

isis = ipi_distribution(data)
ax_isi.hist(isis, bins=20, density = True, label="$\mu(I) = {:.2f}$ \n $C_V(I) = {:.2f}$".format(np.mean(isis), np.std(isis)/np.mean(isis)))
ax_isi.plot(ts, ps, c="C3")
#ax_isi.set_yscale("log")
ax_isi.legend()

amps = amp_distribution(data)
ax_amp.hist(amps, 5, density=True)
ax_amp.set_ylim([0, 1])
plt.savefig(home + "/Data/Calcium/Plots/single_cluster_meanI{:.1f}_ncls{:d}_ncha{:d}.pdf".format(meanI, NClosedStates, NCha), transparent=True)
plt.show()