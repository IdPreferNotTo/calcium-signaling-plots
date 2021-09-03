import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
data = np.loadtxt(file_str)

ISIs = []
MEANs = []
STDs = []
n = len(data[0])
for j in range(1, n):
    ts = [x[0] for x in data]
    cas = [x[j] for x in data]

    spiking: bool = False
    t_tmp: float = 0
    spike_times = []
    if j==2:
        for t, ca in zip(ts, cas):
            if t > 500 and t < 3700:
                if ca > 0.35 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.35 and spiking:
                    spiking = False
    else:
        for t, ca in zip(ts, cas):
            if t > 500 and t < 3700:
                if ca > 0.4 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.4 and spiking:
                    spiking = False

    isis = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]
    n = len(isis)
    mean = np.mean(isis)
    MEANs.append(mean)
    std = np.std(isis)
    STDs.append(std)
    cv = std / mean
    cv2 = np.power(cv, 2)
    print(j, n, mean, std, cv)

    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    for ti in spike_times:
        ax1.axvline(ti, ls=":", c="k", lw=1)
    ax1.plot(ts, cas)
    ax1.axhline(0.4, ls=":", c="C7", lw=1)
    #ax1.set_xlim([0, 3700])
    ax1.set_ylim([0.2, 0.7])

    tmax = max(isis)
    ax2.hist(isis, bins=10, density=True)
    ts = np.linspace(1, tmax, 100)
    inverse_gaussian = []
    for t in ts:
        t3 = np.power(t, 3)
        x = np.sqrt(mean/(2*np.pi*cv2*t3))*np.exp(-np.power(t - mean, 2)/(2*mean*cv2*t))
        inverse_gaussian.append(x)
    ax2.plot(ts, inverse_gaussian)
    plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/ISI/HEK2_bapta_{:d}_isi_dist.png".format(j))
    plt.show()
    plt.close()

fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.scatter(MEANs, STDs)
print(MEANs)
print(STDs)
coef = np.polyfit(MEANs, STDs,1)
print(coef)
t_ref = -coef[1]/coef[0]
print(t_ref)
xs = np.linspace(0, 400, 100)
ys = coef[1] + coef[0]*xs
ax.plot(xs, ys, ls=":", lw=1, c="k")
ax.set_ylim([0, 100])
ax.set_xlim([0, 400])

plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/HEK2_bapta_mean_vs_std.png")
plt.show()
plt.close()


fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

ISIs_no_ref = []
for isis in ISIs:
    isis_no_ref = []
    for isi in isis:
        isis_no_ref.append(isi - t_ref)
    ISIs_no_ref.append([t/np.mean(isis_no_ref) for t in isis_no_ref])

ISIs_flat = [item for sublist in ISIs_no_ref for item in sublist]

mean = np.mean(ISIs_flat)
std = np.std(ISIs_flat)
var = np.power(std, 2)
cv = std / mean
cv2 = np.power(cv, 2)
print(n, mean, std, cv)
tmax = max(ISIs_flat)
ax.hist(ISIs_flat, bins=25, density=True)
ts = np.linspace(0.1, tmax, 100)
inverse_gaussian = []
gaussian = []
for t in ts:
    t3 = np.power(t, 3)
    x = np.sqrt(mean / (2 * np.pi * cv2 * t3)) * np.exp(-np.power(t - mean, 2) / (2 * mean * cv2 * t))
    y = np.sqrt(1./(2*np.pi*var))*np.exp(-np.power(t - mean, 2)/(2*var))
    inverse_gaussian.append(x)
    gaussian.append(y)
ax.plot(ts, inverse_gaussian)
ax.plot(ts, gaussian)

plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/HEK2_bapta_norm_isi_dist.png")
plt.show()
