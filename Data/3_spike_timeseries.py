import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/hek spikes/hek_ca_spikes.dat"
data = np.loadtxt(file_str)

n = len(data[0])
print(n)
for j in range(1, n):
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx-i), axis=0)

    n_avr = 20
    stat_cas = []
    times = [x[0]/1000 for x in data]
    cas = [x[j] for x in data]
    puff_start =[]
    puff_end = []
    avr_cas = []
    puff = False

    for idx, (t, y) in enumerate(zip(times, cas)):
        if idx < n_avr or idx >= (len(cas)-n_avr):
            continue
        mov_avr = np.mean(cas[idx-n_avr: idx])
        mov_std = np.std(cas[idx-n_avr: idx])
        avr_cas.append(mov_avr)
        stat_cas.append(y-mov_avr)
        if y > 1.2*mov_avr and puff == False:
            puff_start.append(t)
            puff = True
        if y < 1.2*mov_avr and puff == True:
            puff = False
            puff_end.append(t)

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1)
    ax2 = fig.add_subplot(gs[1:3, 0])
    ax1 = fig.add_subplot(gs[0, 0], sharex=ax2)
    axis = [ax1, ax2]

    #cas_unbias = [ca - avr_ca for ca, avr_ca in zip(cas, avr_cas)]
    for ptime in puff_start:
        ax1.axvline(ptime)
    ax2.plot(times[n_avr:-n_avr], cas[n_avr:-n_avr])
    ax2.plot(times[n_avr:-n_avr], avr_cas, c="k")

    ax2.plot(times[n_avr:-n_avr], [1.2*x for x in avr_cas], c="C7")
    #ax2.set_ylim([80, 200])
    #ax2.set_xlim([0,50])
    plt.savefig(home + "/Desktop/Ca data/hek spikes/plots/hek spikes idx{:d}.pdf".format(j), transparent=True)
    plt.show()
    plt.close()
