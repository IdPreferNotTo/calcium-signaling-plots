import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc

def transient_func(t, T0, T8, tau):
    return T0*np.exp(-t/tau) + T8*(1 - np.exp(-t/tau))


if __name__ == "__main__":



    mean_intervals_renew = []
    firing_rates_renew = []
    mean_intervals_refrac = []
    firing_rates_refrac = []
    p_effs = []
    tau_effs = []
    tau_effs_theory = []
    tau_cas = []

    iota_cas = []
    iota_effs = []
    iota_eff_stds = []
    iota_effs_theory = []

    ps = [(i/100)*0.2 for i in range(100)]
    tau = 2.0
    delta  = 0.05
    home = os.path.expanduser("~")

    taua = 100
    folder_ref = "/Data/calcium_spikes_markov/Data_adap/"
    for p in ps:
        file = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{delta:.2e}_tau{tau:.2e}_j{p:.2e}_N10_0.dat"
        intervals = np.loadtxt(home + folder_ref + file)
        mean_interval = np.mean(intervals[100:])
        mean_intervals_refrac.append(mean_interval)
        firing_rates_refrac.append(1./mean_interval)
        if len(intervals) != 0:
            #times = [sum(intervals[:i]) for i in range(len(intervals))]
            popt, pcov = curve_fit(transient_func, np.arange(len(intervals)), intervals, p0=(100, 150, 10))
            p_effs.append(p)
            iota_cas.append(taua/mean_interval)
            iota_effs.append(popt[2])
            iota_eff_stds.append(np.sqrt(pcov[2, 2]))
            iota_effs_theory.append(-1/np.log((1-delta)*np.exp(-mean_interval/taua)))

    folder_renew ="/Data/calcium_spikes_markov/Data_no_adap/"
    for p in ps:
        file = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{p:.2e}_N10_0.dat"
        intervals = np.loadtxt(home + folder_renew + file)
        mean_interval = np.mean(intervals)
        mean_intervals_renew.append(mean_interval)
        firing_rates_renew.append(1. / mean_interval)

    intervals_ax2 = np.loadtxt(home + folder_ref + f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{delta:.2e}_tau{tau:.2e}_j1.60e-01_N10_0.dat")
    interval_ax2_0 = np.mean(np.loadtxt(home + folder_renew + f"spike_times_markov_ip1.00_tau{tau:.2e}_j1.60e-01_N10_0.dat"))
    interval_ax2_inf = np.mean(intervals_ax2[100:])

    times = [sum(intervals_ax2[:i]) for i in range(len(intervals_ax2))]
    popt, pcov = curve_fit(transient_func, times, intervals_ax2, p0=(100, 150, 10))
    ts = np.linspace(0, times[30])
    ys = []
    for t in ts:
        y = interval_ax2_0*np.exp(-t/taua) + interval_ax2_inf*(1 - np.exp(-t/taua))
        ys.append(y)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    ax1.set_xlabel("$t$ / s")
    ax1.set_ylabel("$c_i$")

    ax2.set_xlabel("$t$ / s")
    ax2.set_ylabel("$r$")
    ax2.set_ylim([0, 1])
    ax2.plot([sum(intervals_ax2[:i]) for i in range(30)], [1/I for I in intervals_ax2[:30]])
    ax2.plot(ts, [1/y for y in ys])

    ax3.set_xlabel("$p$")
    ax3.set_ylabel("$r$ / s$^{-1}$")
    ax3.plot(ps, firing_rates_refrac)
    ax3.plot(ps, firing_rates_renew, ls=":")

    ax4.set_xlabel("$p$")
    ax4.set_ylabel(r"$\iota$")
    ax4.plot(p_effs, iota_cas)
    ax4.plot(p_effs, iota_effs)
    ax4.fill_between(p_effs, [mean-std for mean, std in zip(iota_effs, iota_eff_stds)], [mean+std for mean, std in zip(iota_effs, iota_eff_stds)], color="C1", alpha=0.3)

    ax4.plot(p_effs, iota_effs_theory)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig4.png", transparent=True)
    plt.show()