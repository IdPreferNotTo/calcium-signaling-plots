import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

import styles as st

home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
data = np.loadtxt(file_str)

n = len(data[0])
for j in range(1, n):
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    print(idxs)
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx-i), axis=0)

    n_avr = 50
    stat_cas = []
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
            if y > 0.35 and not spiking:
                spike_times.append(t)
                spiking = True
            if y < 0.35 and spiking:
                spiking = False

    st.set_default_plot_style()
    fig = plt.figure(tight_layout = True, figsize=(4, 6/2))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    ax1.set_xlim([0, 9000])
    ax1.axhline(0.35, ls = ":", c="C7")
    st.remove_top_right_axis([ax1])

    ISIs = []
    for t1, t2 in zip(spike_times[:-1], spike_times[1:]):
        ISIs.append(t2-t1)
    isi_mean = np.mean(ISIs)
    isi_var =  np.var(ISIs)
    Cv = np.sqrt(isi_var)/isi_mean

    #for spike_time in spike_times:
    #    ax1.axvline(spike_time)
    ax1.plot(t_plot, ca_plot, c=st.colors[4])
    ax1.set_ylabel("Ratio (340/380) [a.u.]")
    ax1.set_xlabel("t [s]")
    plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/Traces/HEK2_bapta_{:d}.png".format(j))
    plt.show()
    plt.close()
