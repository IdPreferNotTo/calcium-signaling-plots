import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

def pre_spike_average(cass):
    n_max = max([len(cas) for cas in cass])
    ca_mean = []
    for n in range(n_max):
        cas_trans= []
        for cas in cass:
            if len(cas) > n:
                cas_trans.append(cas[n])
        ca_mean.append(np.mean(cas_trans))
    return ca_mean


home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/hek spikes/hek_ca_spikes.dat"
data = np.loadtxt(file_str)
data = data[0:-20]

n = 17
row = [x[n] for x in data]
n_avr = 20
stat_cas = []
times = [x[0] for x in data]
cas = [x[n] for x in data]
puff_start = []
puff_start_idx = []
puff_stop = []
puff_stop_idx = []
avr_cas = []
puff = False

for idx, (t, y) in enumerate(zip(times, cas)):
    if idx < n_avr or idx >= (len(cas) - n_avr):
        continue
    mov_avr = np.mean(cas[idx - n_avr: idx])
    mov_std = np.std(cas[idx - n_avr: idx])
    avr_cas.append(mov_avr)
    stat_cas.append(y - mov_avr)
    if y > 1.2 * mov_avr and puff == False:
        puff_start.append(t)
        puff_start_idx.append(idx)
        puff = True
    if y < 1.2 * mov_avr and puff == True:
        puff_stop.append(t)
        puff_stop_idx.append(idx)
        puff = False

pre_spikes = []
pre_spike = []
n = 1
for t, ca in zip(times[::-1], cas[::-1]):
    if t < puff_start[-n]:
        pre_spike.append(ca)
    if n < len(puff_stop_idx):
        if t < puff_stop[-(n+1)]:
            n += 1
            pre_spikes.append(list(pre_spike))
            pre_spike = []
    else:
        if t == 0:
            n += 1
            pre_spikes.append(list(pre_spike))
            pre_spike = []

fig = plt.figure()
gs = gridspec.GridSpec(3, 1)
ax2 = fig.add_subplot(gs[1:3, 0])
ax1 = fig.add_subplot(gs[0, 0], sharex=ax2)
axis = [ax1, ax2]
axin = inset_axes(ax2, width="40%", height="40%", loc=1)

for ps in pre_spikes:
    axin.plot(ps[:-3], c="C7")


# cas_unbias = [ca - avr_ca for ca, avr_ca in zip(cas, avr_cas)]
for pstart, pstop in zip(puff_start, puff_stop):
    ax1.axvline(pstart, c="C4")
    ax1.axvline(pstop, c="C3")
ax2.plot(times[n_avr:-n_avr], cas[n_avr:-n_avr])
ax2.plot(times[n_avr:-n_avr], avr_cas, c="k")

ax2.plot(times[n_avr:-n_avr], [1.2 * x for x in avr_cas], c="C7")
# ax2.set_ylim([80, 200])
# ax2.set_xlim([0,50])
plt.show()
plt.close()
