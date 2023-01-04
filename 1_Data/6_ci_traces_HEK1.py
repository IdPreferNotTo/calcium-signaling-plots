import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

import styles as st

home = os.path.expanduser("~")
file_str = home + "/Data/calcium_spikes_experimental/Spikes/HEK/HEK1_ratio.dat"
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

    st.set_default_plot_style()
    fig = plt.figure(tight_layout = True, figsize=(4, 6/2))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax1])

    #cas_unbias = [ca - avr_ca for ca, avr_ca in zip(cas, avr_cas)]
    ax1.plot(times, cas, lw=1, color=st.colors[4])
    ax1.set_ylim([0,1])
    ax1.set_ylabel("Ratio (340/380) [a.u.]")
    ax1.set_xlabel("t [s]")
    plt.show()
    plt.close()
