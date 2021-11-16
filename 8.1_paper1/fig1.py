import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc

def transient_func2(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))

if __name__ == "__main__":
    home = os.path.expanduser("~")
    file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
    data = np.loadtxt(file_str)
    n = len(data[0])
    j = 12
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    print(idxs)
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx - i), axis=0)
    caT = 0.4
    n_avr = 50
    stat_cas = []
    times = [x[0] for x in data]
    cas = [x[j] for x in data]
    spiking = False
    spike_times = []
    spike_end_times = []
    t_plot = []
    ca_plot = []
    for idx, (t, y) in enumerate(zip(times, cas)):
        if t < 3500:
            if y > caT and not spiking:
                spike_times.append(t)
                spiking = True
            if y < caT and spiking:
                spike_end_times.append(t)
                ca_tmp = []
                t_tmp = []
                spiking = False
            t_plot.append(t)
            ca_plot.append(y)

    ISIs = []
    for t1, t2 in zip(spike_times[:-1], spike_times[1:]):
        ISIs.append(t2 - t1)
    isi_mean = np.mean(ISIs)
    isi_var = np.var(ISIs)
    Cv = np.sqrt(isi_var) / isi_mean

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])

    ax1.set_xlabel("$t$ / s")
    ax1.set_xlim([0, 1500])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])
    ax1.plot(t_plot, ca_plot, c=st.colors[4])
    ax1.set_ylabel("Ratio (340/380) / a.u.")

    nr_ISIs = len(ISIs)
    index_ISIs = np.arange(nr_ISIs)
    popt, pcov = curve_fit(transient_func2, index_ISIs, ISIs, p0=(100, 150, 2))
    print(popt)
    i_ISI_fit = np.linspace(0, nr_ISIs)
    ISI_fit = popt[0]*np.exp(-i_ISI_fit/popt[2]) + popt[1]*(1 - np.exp(-i_ISI_fit/popt[2]))

    ax2.set_xlim([0, nr_ISIs])
    ax2.set_xlabel("$i$")
    ax2.set_ylim([0, 200])
    ax2.set_ylabel("$T_i$ / s")
    ax2.scatter(index_ISIs, ISIs, fc="w", ec=st.colors[4], s=20, zorder=3)
    ax2.plot(i_ISI_fit, ISI_fit, lw=1, c="k", zorder=2)
    ax2.axhline(popt[1], ls=":", lw=1, c="k")
    ax2.text(nr_ISIs+2, popt[1], "$T_0$", va="center", ha="center")
    ax2.fill_between([0, 10], [0, 0], [200, 200], color="C7", alpha=0.5, zorder=1)

    ax3.set_xlabel("$t$ / s")
    ax3.set_ylabel("$c_i$")


    ax4.set_xlabel("$i$")
    ax4.set_ylabel("$T_i$ / s")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig1_ca_integrate_and_fire.pdf",
                transparent=True)
    plt.show()
    plt.close()
