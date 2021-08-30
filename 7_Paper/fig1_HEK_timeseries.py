import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import functions as fc

home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
data = np.loadtxt(file_str)

ISIs = []
ISIs_no_transient = []
MEANs = []
STDs = []
n = len(data[0])

j = 12
ts = [x[0] for x in data]
cas = [x[j] for x in data]

spiking: bool = False
t_tmp: float = 0
spike_times = []
for t, ca in zip(ts, cas):
    if t < 3500:
        if ca > 0.4 and not spiking:
            spike_times.append(t)
            spiking = True
        if ca < 0.4 and spiking:
            spiking = False


ISIs = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]
ISIs_no_transient = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1]) if t1 > 500]
n = len(ISIs)
mean = np.mean(ISIs)
std = np.std(ISIs)

cv = std / mean
cv2 = np.power(cv, 2)
print(j, n, mean, std, cv)

fc.set_default_plot_style()
fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
gs = gridspec.GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0,:])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
fc.remove_top_right_axis([ax1, ax2, ax3])
for ti in spike_times:
    ax1.axvline(ti, ymin=0.0, ymax=0.8, ls=":", c="k", lw=1)
ax1.plot(ts, cas, lw=1)
ax1.set_xlim([0, 1000])
ax1.set_ylim([0.2, 1.0])
ax1.set_xlabel("$t$ / [s]")
ax1.set_ylabel("Ratio (340/380)")

for i, (ti2, ti1) in enumerate(zip(spike_times[1:], spike_times[:-1])):
    if ti2 < 1000:
        ax1.text((ti1 + ti2)/2, 0.9, f"$T_{i}$", ha="center", va="center", clip_on=False)

        ax1.arrow(ti1, 0.8, ti2-ti1, 0, fc="k", length_includes_head=True, head_width=0.05, head_length=25.0,
             lw=0.5,
             clip_on=False)
        ax1.arrow(ti2, 0.8, ti1 -ti2, 0, fc="k", length_includes_head=True, head_width=0.05, head_length=25.0,
             lw=0.5,
             clip_on=False)

T = 0
for n, t in enumerate(spike_times):
    if t > 500:
        break

ax2.set_xlim(0, len(ISIs))
ax2.set_ylim([0, 1.1*max(ISIs)])
ax2.set_xlabel("$k$")
ax2.set_ylabel("$T_k$ / [s]")
ax2.set_xticks(range(0, len(ISIs), 10))
ks = np.arange(len(ISIs))
ax2.axvline(n, ls=":", c="k", lw=1)
ax2.scatter(ks, ISIs, s=20, fc="w", ec="C0")
ax2.axvspan(0, n, facecolor="C7", alpha=0.5, zorder=0)

ax3.set_xlim(0.8*min(ISIs_no_transient), 1.1*max(ISIs_no_transient))
ax3.hist(ISIs_no_transient, bins=range(int(min(ISIs_no_transient)), int(max(ISIs_no_transient)) + 5, 5))
#ax3.set_xlim([100, 200])
ax3.set_xlabel("$T_k$ / [s]")
ax3.set_ylabel("$P(T_k)$")

ax1.text(0.5, 1.0, "A B C", size=12, transform=ax1.transAxes)
plt.savefig(home + f"/Paper/RamLin21_Ca/Plots/fig1_HEK_timeseries_j{j}.pdf")
plt.show()
plt.close()


