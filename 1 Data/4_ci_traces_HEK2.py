import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

import styles as st

home = os.path.expanduser("~")
file_str = home + "/Data/calcium/experimental/Spikes/HEK/HEK2_bapta_ratio.dat"
data = np.loadtxt(file_str)

thresholds = [0.4, 0.35, 0.4, 0.35, 0.4, 0.35, 0.35, 0.35, 0.35, 0.4, 0.45, 0.35, 0.35, 0.45, 0.4, 0.35, 0.4, 0.45,
              0.35, 0.35, 0.5, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.4, 0.45, 0.45, 0.45,
              0.35, 0.4, 0.4]

n = len(data[0])
for j in range(1, n):
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    print(idxs)
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx-i), axis=0)

    times = [x[0] for x in data]
    cas = [x[j] for x in data]

    spiking = False
    spike_times = []
    t_plot = []
    ca_plot= []
    for idx, (t, y) in enumerate(zip(times, cas)):
        if t < 10_000:
            t_plot.append(t)
            ca_plot.append(y)
            if y > thresholds[j-1] and not spiking:
                spike_times.append(t)
                spiking = True
            if y < thresholds[j-1] and spiking:
                spiking = False

    st.set_default_plot_style()
    fig = plt.figure(tight_layout = True, figsize=(4, 2))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax1])

    ISIs = []
    for t1, t2 in zip(spike_times[:-1], spike_times[1:]):
        ISIs.append(t2-t1)
    isi_mean = np.mean(ISIs)
    isi_var =  np.var(ISIs)
    Cv = np.sqrt(isi_var)/isi_mean

    #for spike_time in spike_times:
    #    ax1.axvline(spike_time)
    ax1.set_xlim([0, 4000])
    ax1.set_ylim([0.2, 1])
    #ax1.axhline(thresholds[j-1], ls = "--", c="C7")
    #for t in spike_times:
        #ax1.axvline(t, lw=1, c="C7")
    ax1.plot(t_plot, ca_plot, c=st.colors[7], lw=1)
    ax1.set_ylabel("Ratio (340/380)")
    ax1.set_xlabel("$t$ / s")
    plt.savefig(
        home + "/Data/calcium/experimental/Spikes/HEK/Plots2/Traces/HEK2_bapta_{:d}.png".format(j), dpi=300, transparent=True)

    plt.show()
    plt.close()
