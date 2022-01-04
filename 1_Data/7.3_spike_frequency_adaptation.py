import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from scipy.optimize import curve_fit

import functions as fc
import styles as st

def transient_func(t, T0, T8, tau):
    return T0*np.exp(-t/tau) + T8*(1 - np.exp(-t/tau))

home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
data = np.loadtxt(file_str)

CVs = []
T0s = []
Tinfs = []
TAUs = []
n = len(data[0])
for j in range(1, n):
    ts = [x[0] for x in data]
    cas = [x[j] for x in data]

    spiking: bool = False
    t_tmp: float = 0
    spike_times = []
    if j==2:
        for t, ca in zip(ts, cas):
            if t < 3700:
                if ca > 0.35 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.35 and spiking:
                    spiking = False
    else:
        for t, ca in zip(ts, cas):
            if t < 3700:
                if ca > 0.4 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.4 and spiking:
                    spiking = False
    print(j)
    ISIs = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]
    idx_ISIs = [i for i in range(len(ISIs))]
    t_ISIs = []
    T = 0
    for isi in ISIs:
        T += isi
        t_ISIs.append(T)
    #popt, pcov = curve_fit(transient_func, t_ISIs, ISIs, p0=(10, 100, 200))
    popt, pcov = curve_fit(transient_func, idx_ISIs, ISIs, p0=(10, 100, 10))
    T0s.append(popt[0])
    Tinfs.append(popt[1])
    TAUs.append(popt[2])
    CVs.append(np.std(ISIs)/np.mean(ISIs))

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax1])

    print(popt,  np.std(ISIs)/np.mean(ISIs))
    #xts = np.linspace(0, 4000)
    #yTs = popt[0]*np.exp(-xts/popt[2]) + popt[1]*(1 - np.exp(-xts/popt[2]))
    Tis = [popt[0] * np.exp(-idx/popt[2]) + popt[1] * (1 - np.exp(-idx/popt[2])) for idx in idx_ISIs]
    ax1.plot(idx_ISIs, Tis, c="k", lw=1, label=f"$f_0$ = {popt[0]:.1e}" + "\n" + f"$f_\infty$ = {popt[1]:.1e}" + "\n" + rf"$\tau_{{eff}}$ = {popt[2]:.1f}")
    ax1.set_ylim([0, 1.1*max(ISIs)])
    ax1.set_xlabel("$i$")
    ax1.set_ylabel(r"$T_i$ / s")
    ax1.axvline(10, ls=":", c="k", lw=1)
    ax1.scatter(idx_ISIs, ISIs, s=20, fc="w", ec="C0")
    ax1.axvspan(0, 10, facecolor="C7", alpha=0.5, zorder=0)
    legend = ax1.legend(loc=4, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 10})
    legend.get_frame().set_linewidth(0.5)
    plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/Adaptation/HEK2_bapta_{:d}_adaptation.png".format(j))
    plt.show()