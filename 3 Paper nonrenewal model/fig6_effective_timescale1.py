import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

def exponential_cer(t, cer0, cer8, tau):
    return cer8 + (cer0 - cer8) * np.exp(-t / tau)


def turn_2d_isis_into_spike_times(ISIss):
    spike_times = np.zeros(ISIss.shape)
    for i, ISIs in enumerate(ISIss):
        t = 0
        for ii, I in enumerate(ISIs):
            t += I
            spike_times[i, ii] = t
    return spike_times


if __name__ == "__main__":
    # Set Plot style
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.25 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(3, 1, height_ratios= [0.2, 0.4, 0.4], hspace=0.)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])
    st.remove_everything([ax0])
    st.remove_top_right_axis([ax1, ax2])

    ax1.set_ylabel(r"$r(t)$ / s\textsuperscript{-1}")
    ax1.set_xticklabels([])
    ax1.set_xlim([-600, 2100])

    ax2.set_ylim([0.875, 1.05])
    ax2.set_xlim([-600, 2100])
    ax2.set_ylabel(r"$c_{\rm er}(t)$")
    ax2.set_xlabel(r"$t$ / s")

    # Parameters
    tau = 5
    p = 0.015
    eps_er = 0.03
    tau_er = 300

    # Plot Stimulation
    ax0.set_ylim([-0.5, 1.5])
    ax0.plot([-500, 0, 0, 1000, 1000, 2000], [0, 0, 1, 1, 0, 0], c="k")
    ax0.text(-250, 0.07, "Stim.\ off", ha="center", va='bottom')
    ax0.text(500, 1.07, "Stim.\ on", ha="center", va='bottom')
    ax0.text(1500, 0.07, "Stim.\ off", ha="center", va='bottom')

    # Plot spike times by vertical lines
    data = df.load_traces_markov(tau, p, True, tau_er, eps_er)
    ts, cs, js, cers = np.transpose(data)
    # Account for off_signal in data from t=0 to t=100

    ax2.plot(ts[1_000:11_000]-100, cers[1_000:11_000], lw=1.15, c="k", zorder=1)
    ts_spike = []
    for t, c, j in zip(ts, cs, js):
        if j == 69:
            ts_spike.append(t)
    for t_spike in ts_spike[:10]:
        ax0.plot([t_spike-100, t_spike-100], [0.2, 0.8], lw=1.15, c="k")

    # Plot r(t)
    data_isi_2d = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er)
    spike_times_2d = turn_2d_isis_into_spike_times(data_isi_2d)
    spike_times_1d = [item for sublist in spike_times_2d for item in sublist]
    spike_times_sorted = sorted(spike_times_1d)
    ts = np.linspace(0, 1000, 100)
    n_over_t = np.zeros(100)
    for t1 in spike_times_sorted:
        i = int(t1/10)
        if i >= 100:
            continue
        else:
            n_over_t[i] += 1
    r_over_t = n_over_t / (10. * 1000)

    ax1.plot([-500, 0], [0, 0], zorder=2, lw=1.15, c="C7")
    ax1.plot(ts, r_over_t, zorder=2, lw=1.15, c="C7")
    ax1.plot([1000, 1000, 2000], [r_over_t[-1], 0, 0], lw=1.15, c="C7")
    ax1.axvline(0, ls=":", c="C7")
    ax1.axvline(1000, ls=":", c="C7")
    formatter = mticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-1, 2))
    ax1.yaxis.set_major_formatter(formatter)

    # Plot c_er(t)
    file = f"transient_adaptation_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
    data = np.loadtxt("/home/lukas/Data/calcium/markov/adaptive/transient/" + file,
                      usecols=np.arange(0, 1000))
    mean_cers_on = np.mean(data, axis=0)
    ts = np.linspace(0, 1000, 1000)
    popt_on, pcov = curve_fit(fc.exponential_cer, ts, mean_cers_on, p0=(0.9, 100))
    mean_cers_on_fit = popt_on[0] + (1 - popt_on[0]) * np.exp(-ts / popt_on[1])

    mean_cer_off1 = mean_cers_on_fit[-1]
    mean_cer_off = mean_cer_off1
    mean_cers_off = []
    dt = ts[1] - ts[0]
    for t in ts:
        mean_cers_off.append(mean_cer_off)
        dcer = -(mean_cer_off - 1.) / tau_er
        mean_cer_off += dcer*dt
    popt_off, pcov = curve_fit(exponential_cer, ts, mean_cers_off, p0=(0.9, 1.0, 100))
    mean_cers_off_fit = popt_off[1] + (popt_off[0] - popt_off[1]) * np.exp(-ts / popt_off[2])

    #ax2.text(-200, 0.9, "Fast depletion \n with" + r"$\tau_{\rm eff}$", ha="center", va='bottom')
    ax2.annotate("Fast depletion \n with " + r"$\tau_{\rm eff}$", xy=(50, 0.965),
                 xytext=(-550, 0.9), va="center", ha="left",
                 arrowprops=dict( arrowstyle="->" ))
    ax2.annotate("Slow replenishment \n with " + r"$\tau_{\rm er}$", xy=(1250, 0.965),
                 xytext=(1200, 0.9), va="center", ha="left",
                 arrowprops=dict( arrowstyle="->" ))

    ax2.axvline(0, ls=":", lw=1.15, c="C7")
    ax2.axvline(1000, ls=":", lw=1.15, c="C7")

    ts = np.insert(ts, 0, -500)
    mean_cers_on = np.insert(mean_cers_on, 0, 1)
    mean_cers_off = np.insert(mean_cers_off, 0, mean_cer_off1)
    ax2.plot(ts, mean_cers_on, color="C7")
    ax2.plot(ts + 1000, mean_cers_off, color="C7")

    cer_fix = fc.self_consistent_cer_infty(tau, p, tau_er, eps_er)
    tau1, _ = fc.calculate_tau_1(tau, p, tau_er, eps_er)
    tau2, _ = fc.calculate_tau_2(tau, p, tau_er, eps_er)
    ts = np.linspace(0, 1000, 1000)
    cers1 = cer_fix + (1 - cer_fix) *  np.exp(-ts/tau1)

    cers2 = cer_fix + (1 - cer_fix) *  np.exp(-ts/tau2) #cer_fix * ((1 + cer_fix) +(1 - cer_fix)*np.exp(-ts/tau2))/((1 + cer_fix) - (1 - cer_fix)*np.exp(-ts/tau2))
    #ax2.plot(ts, cers1, c=st.colors[7], ls=":")
    ax2.plot(ts, cers2, c=st.colors[7], ls=":") #(0, (3, 1)))

    data_rate_cer = np.loadtxt(home + f"/Data/calcium/theory/cer_r0_fpe_tau{tau:.2e}_j{p:.2e}.dat")
    cer_fix, r_fix = np.transpose(data_rate_cer)
    r1s = []
    for cer in cers1:
        r = fc.linear_interpolate(cer, cer_fix, r_fix)
        r1s.append(r)
    r2s = []
    for cer in cers2:
        r = fc.linear_interpolate(cer, cer_fix, r_fix)
        r2s.append(r)


    #ax1.plot(ts, r1s, c=st.colors[7], ls=":")
    ax1.plot(ts, r2s, c=st.colors[7], ls=":")#(0, (3, 1)))
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig6.pdf", dpi=300, transparent=True)
    plt.show()

