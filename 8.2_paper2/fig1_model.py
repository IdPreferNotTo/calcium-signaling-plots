import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.optimize import curve_fit

import os

import styles as st
import default_parameters as df
import functions as fc

def get_firing_rate(x, r_over_cer):
    cers, r0s = np.transpose(r_over_cer)
    idx = 100
    for i, cer in enumerate(cers):
        if cer >= x:
            idx = i
            break
    return r0s[idx]


def turn_2d_isis_into_spike_times(ISIss):
    spike_times = np.zeros((6000, 25))
    for i, ISIs in enumerate(ISIss):
        t = 0
        for ii, I in enumerate(ISIs):
            t += I
            spike_times[i, ii] = t
    return spike_times


if __name__ == "__main__":
    # Set up plot style
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.66 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(nrows=2, ncols=2)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    axis = [ax1, ax2, ax3, ax4]
    st.remove_top_right_axis(axis)
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"C", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"D", fontsize=11, transform=ax4.transAxes, va='top')

    ca_r = 0.2
    ca_t = 0.5
    ax1.set_xlabel("$t$ / s")
    ax1.set_ylabel(r"$c_{\rm i}(t)$")
    ax1.set_xlim([-100, 600])
    ax1.set_yticks([ca_r, ca_t])
    ax1.set_ylim([0.1, 1.1])
    ax1.set_yticklabels(["$c_R$", "$c_T$"])
    ax1.axhline(ca_t, ls=":", lw=1, c="k")

    ax2.set_xlabel("$t$ / s")
    ax2.set_ylabel(r"$c_{\rm er}(t)$")
    ax2.set_ylim([0.85, 1.05])
    ax2.set_xlim([-100, 600])

    ax3.set_xlabel("$i$")
    ax3.set_ylabel("$T_i$ / s")

    ax4.set_ylabel(r"$p_{\rm ISI}(t)$")
    ax4.set_xlabel(r"$t / \langle T \rangle$")

    tau = 5.0
    p = 0.015
    eps_er = 0.03
    tau_er = 300

    # Plot traces
    home = os.path.expanduser("~")
    data_traces = df.load_traces_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
    ts, cis, jps, cers = np.transpose(data_traces)
    ts_plot = []
    cis_plot = []
    cers_plot = []
    n_spikes = 0
    spike_times = []
    cr_spike_times = []

    for t, ci, jp, cer in zip(ts, cis, jps, cers):
        t = t - 100
        if t >= 1_000:
            break
        ts_plot.append(t)
        cis_plot.append(ci)
        cers_plot.append(cer)

        if jp == 69:
            spike_times.append(t)
            ts_plot.append(t)
            cis_plot.append(0.5 + cer / 2)
            cers_plot.append(cer)
            ts_plot.append(t)
            cis_plot.append(0.2)
            cers_plot.append(cer)
            cr_spike_times.append(0.2*cer)
            n_spikes += 1

    ax1.plot(spike_times, cr_spike_times, ls=":", lw=1, c="k")
    ax1.plot(ts_plot, cis_plot, lw=1, color=st.colors[0])
    ax2.plot(ts_plot, cers_plot, lw=1, color=st.colors[0], label = r"$c_{\rm er}(t)$")

    #ax2.text(300, 1.05, r"$\{ c_{{\rm er}, i}\}$", fontsize=8, ha="center", va='center', bbox=dict(boxstyle="round", fc="w", ec="w"))
    #for i in range(7):
    #    p1 = patches.FancyArrowPatch((300, 1.05), (spike_times[i] , ceris[i]), lw=0.5, arrowstyle='->',
    #                                 mutation_scale=5)
    #    ax2.add_patch(p1)

    # Plot time-dependent rate and Ti
    data_isi_trans = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er)
    spike_times_trans = turn_2d_isis_into_spike_times(data_isi_trans)
    mean_isis = np.mean(data_isi_trans, axis=0)

    data_isi = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er)
    rows, cols = data_isi.shape
    idx_max = 20
    idxs = np.arange(idx_max)
    means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs]  # calculate the mean column wise

    for row in np.arange(1):
        ax3.scatter(idxs, data_isi[row, :idx_max], fc="w", ec=st.colors[1], alpha=0.75, s=20, zorder=1)
    ax3.scatter(idxs, means_Tidx, fc="w", ec="C7", s=20, zorder=2)
    popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=(100, 150, 2))
    ax3.axvspan(0, popt[2], alpha=0.3, color="C7")

    # print(popt)
    ISI_fit = popt[0] * np.exp(-idxs / popt[2]) + popt[1] * (1. - np.exp(-idxs / popt[2]))
    ax3.plot(idxs, ISI_fit, c="k", zorder=3)
    ax3.axhline(popt[0], ls=":", lw=1, c="k")
    ax3.axhline(popt[1], ls=":", lw=1, c="k")
    ax3.set_xlim([0, len(idxs) + 5])
    p1 = patches.FancyArrowPatch((len(idxs) + 1, popt[0]), (len(idxs) + 1, popt[1]), arrowstyle='<->',
                                 mutation_scale=10)
    ax3.add_patch(p1)
    ax3.text(len(idxs) + 2, (popt[0] + popt[1]) / 2, r"$\Delta T$", fontsize=8, ha="left", va='center')

    p1 = patches.FancyArrowPatch((0 - 0.5, popt[1] * 1.1), (popt[2] + 1, popt[1] * 1.1), arrowstyle='<->', mutation_scale=10)
    ax3.add_patch(p1)
    ax3.text(popt[2] / 2, popt[1] * 1.2, r"$n_{\rm tr}$", fontsize=8, ha="center", va='bottom')
    file = home + f"/Data/calcium/markov/adaptive/sequence/long_spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
    data_isi = np.loadtxt(file)
    mean = np.mean(data_isi[20:])
    eps = eps_er / (1 - eps_er / 2)
    cer_fix = 1 / (1 + eps * tau_er / mean)
    ax2.text(0, cer_fix * 1.01, r"$\langle c_{\rm er}^* \rangle$", fontsize=8, ha="center", va='bottom')
    ax2.axhline(cer_fix, ls=":", lw=1, c="k")
    std = np.std(data_isi)
    cv = std/mean
    cv2 = cv**2
    print(mean, cv)
    ts = np.linspace(0, 1+4*cv, 501)
    p_inv_gaus = fc.inverse_gaussian_dist(ts, 1, cv2)
    l1, = ax4.plot(ts, p_inv_gaus, lw=1, c="k", label="Inv.\ Gaus.")
    ax4.hist(data_isi/mean, bins=25, color=st.colors[1], alpha=0.55, density=True)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig1.pdf")
    plt.show()