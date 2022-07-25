import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc

def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))

if __name__ == "__main__":
    home = os.path.expanduser("~")
    file_str = home + "/Data/calcium_spikes_experimental/Spikes/HEK/HEK2_bapta_ratio.dat"
    data = np.loadtxt(file_str)
    n = len(data[0])
    j = 13
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
    for t1, t2 in zip(spike_times[1:-1], spike_times[2:]):
        ISIs.append(t2 - t1)

    Mean = np.mean(ISIs[10:])
    Cv = np.std(ISIs[10:])/Mean
    print(f"Mean: {Mean:.2f}")
    print(f"CV: {Cv:.2f}")

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 6))
    gs = gridspec.GridSpec(3, 2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    ax5 = fig.add_subplot(gs[2, 0])
    ax6 = fig.add_subplot(gs[2, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4, ax5, ax6])

    ax1.text(0.1, 0.95, "A$_{i}$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.075, 0.95, "A$_{ii}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.1, 0.95, "B$_{i}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.075, 0.95, "B$_{ii}$", fontsize=13, transform=ax4.transAxes, va='top')
    ax5.text(0.1, 0.95, "C$_{i}$", fontsize=13, transform=ax5.transAxes, va='top')
    ax6.text(0.075, 0.95, "C$_{ii}$", fontsize=13, transform=ax6.transAxes, va='top')

    ax1.text(0.6, 1.00, "Data", fontsize=15, transform=ax1.transAxes, ha="center", va='top')
    ax5.text(0.6, 1.00, "Cumulative \n refractoriness", fontsize=15, transform=ax5.transAxes, ha="center", va='top')
    ax3.text(0.6, 1.00, "Renewal", fontsize=15, transform=ax3.transAxes, ha="center", va='top')

    ax1.set_xlabel("$t$ / s")
    ax1.set_xlim([0, 1000])
    ax1.plot(t_plot, ca_plot, c=st.colors[4])
    ax1.set_ylabel("Ratio (340/380) / a.u.")

    nr_ISIs = len(ISIs)
    index_ISIs = np.arange(nr_ISIs)
    popt, pcov = curve_fit(transient_func, index_ISIs, ISIs, p0=(100, 150, 2))
    print(popt)
    i_ISI_fit = np.linspace(0, nr_ISIs)
    ISI_fit = popt[0]*np.exp(-i_ISI_fit/popt[2]) + popt[1]*(1 - np.exp(-i_ISI_fit/popt[2]))

    ax2.set_xlim([0, nr_ISIs])
    ax2.set_xlabel("$i$")
    ax2.set_ylim([0, 1.5*popt[1]])
    ax2.set_ylabel("$T_i$ / s")
    ax2.scatter(index_ISIs, ISIs, fc="w", ec=st.colors[4], s=20, zorder=3)
    ax2.plot(i_ISI_fit, ISI_fit, lw=1, c="k", zorder=2)
    ax2.fill_between([0, 10], [0, 0], [200, 200], color="C7", alpha=0.5, zorder=1)
    ax2.axhline(popt[0], ls=":", lw=1, c="k")
    ax2.axhline(popt[1], ls=":", lw=1, c="k")
    ax2.text(nr_ISIs/2, popt[0]*1.2, "$T_0$", ha="center")
    ax2.text(nr_ISIs/2, popt[1]*1.2, "$T_\infty$", ha="center")

    taua = 829
    ampa = 0.0309
    home = os.path.expanduser("~")
    folder = home + "/Data/calcium_spikes_markov/Data_adap/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    file_spikes = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    data = np.loadtxt(folder + file)
    ts, cas, jpuffs, adaps = np.transpose(data)
    ISIs = np.loadtxt(folder + file_spikes)


    ts_plot = []
    cas_plot = []
    n_spikes = 0
    for t, c, a in zip(ts, cas, adaps):
        ts_plot.append(t + n_spikes * 10)
        cas_plot.append(c)
        if c == 1:
            ts_plot.append(t + n_spikes * 10)
            cas_plot.append(1 + 2*a)
            ts_plot.append(t + n_spikes * 10)
            cas_plot.append(0.33)
            n_spikes += 1
    ax5.axhline(1, lw=1, ls=":", color="C7")
    ax5.set_xlim([0, 1000])
    ax5.set_xlabel("$t$ / s")
    ax5.set_ylim([0.2, 4])
    ax5.set_yticks([0.33, 1])
    ax5.set_yticklabels(["$c_R$", "$c_T$"])
    ax5.set_ylabel("$c_i$ / a.u.")
    ax5.plot(ts_plot, cas_plot, lw=1, color=st.colors[1])

    ax6.set_xlim([0, nr_ISIs])
    ax6.set_xlabel("$i$")
    ax6.set_ylim([0, 1.5 * popt[1]])
    ax6.set_ylabel("$T_i$ / s")
    ax6.scatter(np.arange(0, nr_ISIs), ISIs[0:nr_ISIs], fc="w", ec=st.colors[1], s=20, zorder=3)
    ax6.fill_between([0, 10], [0, 0], [200, 200], color="C7", alpha=0.5, zorder=1)
    index_ISIs = np.arange(len(ISIs))
    popt, pcov = curve_fit(transient_func, index_ISIs, ISIs, p0=(100, 150, 2))
    print(popt)
    ISI_fit = popt[0] * np.exp(-index_ISIs / popt[2]) + popt[1] * (1. - np.exp(-index_ISIs / popt[2]))
    ax6.plot(index_ISIs, ISI_fit, lw=1, c="k")
    ax6.axhline(popt[0], ls=":", lw=1, c="k")
    ax6.axhline(popt[1], ls=":", lw=1, c="k")
    ax6.text(nr_ISIs/2, popt[0]*1.2, "$T_0$", ha="center")
    ax6.text(nr_ISIs/2, popt[1]*1.2, "$T_\infty$", ha="center")

    folder = home + "/Data/calcium_spikes_markov/Data_no_adap/"
    file = "ca_markov_ip1.00_tau1.05e+01_j1.10e-02_N10_0.dat"
    file_spikes = "spike_times_markov_ip1.00_tau1.05e+01_j1.10e-02_N10_0.dat"
    data = np.loadtxt(folder + file)
    ts, cas, jpuffs, adaps = np.transpose(data)
    ISIs = np.loadtxt(folder + file_spikes)

    ts_plot = []
    cas_plot= []
    n_spikes = 0
    for t, c in zip(ts, cas):
        ts_plot.append(t + n_spikes*10)
        cas_plot.append(c)
        if c == 1:
           ts_plot.append(t + n_spikes*10)
           cas_plot.append(3)
           ts_plot.append(t + n_spikes*10)
           cas_plot.append(0.33)
           n_spikes += 1

    ax3.axhline(1, lw=1, ls=":", color="C7")
    ax3.set_xlim([0, 1000])
    ax3.set_xlabel("$t$ / s")
    ax3.set_ylim([0.2, 4])
    ax3.set_yticks([0.33, 1])
    ax3.set_yticklabels(["$c_R$", "$c_T$"])
    ax3.set_ylabel("$c_i$ / a.u.")
    ax3.plot(ts_plot, cas_plot, lw=1, color=st.colors[1])

    ax4.set_xlim([0, nr_ISIs])
    ax4.set_xlabel("$i$")
    ax4.set_ylim([0, 1.5 * popt[1]])
    ax4.set_ylabel("$T_i$ / s")
    ax4.scatter(np.arange(0, nr_ISIs), ISIs[0:nr_ISIs], fc="w", ec=st.colors[1], s=20, zorder=3)
    ax4.axhline(popt[1], ls=":", lw=1, c="k")
    ax4.text(nr_ISIs / 2, popt[1] * 1.2, "$T_0 = T_\infty$", ha="center")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig1.png")
    plt.show()
    plt.close()
