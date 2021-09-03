import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK1_ratio.dat"
data = np.loadtxt(file_str)

n = len(data[0])
for j in range(1, n):
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    print(idxs)
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx-i), axis=0)

    n_avr = 50
    times = [x[0] for x in data]
    cas = [x[j] for x in data]
    spike_times =[]
    spike = False

    for idx, (t, ca) in enumerate(zip(times, cas)):
        if ca > 0.4 and spike == False:
            spike_times.append(t)
            spike = True
        if ca < 0.4:
            spike = False

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1)
    ax2 = fig.add_subplot(gs[1:3, 0])
    ax1 = fig.add_subplot(gs[0, 0], sharex=ax2)
    axis = [ax1, ax2]
    axin = inset_axes(ax2, width="40%", height="40%", loc=1)
    ipis = []
    for t1, t2 in zip(spike_times[:-1], spike_times[1:]):
        ipis.append(t2-t1)
    ipi_mean = np.mean(ipis)
    ipi_var =  np.var(ipis)
    Cv = np.sqrt(ipi_var)/ipi_mean

    #cas_unbias = [ca - avr_ca for ca, avr_ca in zip(cas, avr_cas)]
    for spike_time in spike_times:
        ax1.axvline(spike_time)
    ax2.plot(times, cas)
    ax2.set_ylim([0,1])
    plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots1/HEK1_{:d}.png".format(j))
    plt.show()
    plt.close()
