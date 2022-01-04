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
j = 8
row = [x[j] for x in data]
idxs = [i for (i, x) in enumerate(row) if x > 500]
print(idxs)
for i, idx in enumerate(idxs):
    data = np.delete(data, (idx - i), axis=0)

n_avr = 50
stat_cas = []
times = [x[0] for x in data]
cas = [x[j] for x in data]

caT = 0.32
caR = 0.25

spiking = False
spike_times = []
t_plot = []
ca_plot = []

ca_tmp = []
t_tmp = []
cas_plot = []
ts_plot = []
for idx, (t, y) in enumerate(zip(times, cas)):
    if t < 3000:
        if y > caT and not spiking:
            t_tmp.append(t)
            ca_tmp.append(caT)
            ts_plot.append(t_tmp.copy())
            cas_plot.append(ca_tmp.copy())
            spike_times.append(t)
            spiking = True
        if y < caT and spiking:
            ca_tmp = []
            t_tmp = []
            spiking = False
        t_plot.append(t)
        ca_plot.append(y)
        t_tmp.append(t)
        ca_tmp.append(y)

ISIs = []
for t1, t2 in zip(spike_times[:-1], spike_times[1:]):
    ISIs.append(t2 - t1)
isi_mean = np.mean(ISIs)
isi_var = np.var(ISIs)
Cv = np.sqrt(isi_var) / isi_mean

st.set_default_plot_style()
fig = plt.figure(tight_layout=True, figsize=(4, (3/2)*6 / 2))
gs = gridspec.GridSpec(2, 1)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax1.set_xlim([0, 3000])
ax1.set_ylim([0.2, 0.6])
ax1.axhline(caT, lw=1, c="C7")
ax1.axhline(caR, lw=1, c="C7")
st.remove_top_right_axis([ax1, ax2])
ax1.plot(t_plot, ca_plot, c=st.colors[4])
ax1.set_ylabel("Ratio (340/380) [a.u.]")

ax2.set_xlim([0, 3000])
ax2.set_xlabel("t [s]")
ax2.set_ylabel("$c_i$ [a.u.]")
ax2.set_ylim([0.2, 0.6])
ax2.axhline(caT, lw=1, c="C7")
ax2.axhline(caR, lw=1, c="C7")
ax2.text(3020, caR, "$c_R$", va="center")
ax2.text(3020, caT, "$c_T$", va="center")

for (cas, ts) in zip(cas_plot, ts_plot):
    ax2.plot([ts[2], ts[2]], [cas[2], caT], ls=":", c="C7")
    ax2.plot(ts[2:], cas[2:], c=st.colors[4])
    ax2.arrow(x=ts[-1], y=caT, dx=0, dy=0.6 - caT, color="k", length_includes_head=True, head_width=36., head_length=0.04)


t_start = ts_plot[5][-1]
t_ref = ts_plot[6][2]
t_end = ts_plot[6][-1]

shift = 2.5
ax2.plot([t_start+shift, t_start+shift], [0.36, 0.38],zorder=2, c="C3")
ax2.plot([t_ref-shift, t_ref-shift], [0.36, 0.38],zorder=2, c="C3")
ax2.plot([t_ref+shift, t_ref+shift], [0.36, 0.38],zorder=1, c="C7")
ax2.plot([t_end-shift, t_end-shift], [0.36, 0.38],zorder=1, c="C7")

ax2.plot([t_start, t_ref], [0.37, 0.37], zorder=2, c="C3")
ax2.plot([t_ref, t_end], [0.37, 0.37], zorder=1, c="C7")

plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/Traces/HEK2_bapta_if_dynamics{:d}.png".format(j))
plt.show()
plt.close()
