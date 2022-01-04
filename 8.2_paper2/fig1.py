import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import gridspec
import os

import styles as st


def transient_func(i, T0, T8, tau):
    return T0 * np.exp(-i / tau) + T8 * (1 - np.exp(-i / tau))


if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 3))
    gs = gridspec.GridSpec(1, 2)
    gs1 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0])
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.75)
    ax0 = fig.add_subplot(gs1[0])
    ax1 = fig.add_subplot(gs1[1])
    ax2 = fig.add_subplot(gs1[2])
    ax3 = fig.add_subplot(gs2[0])
    ax4 = fig.add_subplot(gs2[1])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3, ax4])

    taua = 829
    ampa = 0.0309
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    file_spikes = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    data = np.loadtxt(folder + file)
    ts, cas, jpuffs, adaps = np.transpose(data)
    ISIs = np.loadtxt(folder + file_spikes)

    mean_isi = np.mean(ISIs)
    cv_isi = np.std(ISIs) / mean_isi
    print(mean_isi, cv_isi)

    spike_times = []
    cR = 0.33
    cT = 1.0
    t_isi = []
    ca_isi = []
    jpuff_isi = []
    adap_isi = []
    nrISIs = 0
    for t, ca, jpuff, adap, adap_after in zip(ts[:-1], cas[:-1], jpuffs[:-1], adaps[:-1], adaps[1:]):
        t_isi.append(t)
        ca_isi.append(ca)
        jpuff_isi.append(jpuff)
        adap_isi.append(adap)
        if ca == 1 and nrISIs <= 10:
            nrISIs += 1
            spike_times.append(t)

            ax0.plot(t_isi, ca_isi, lw=1, c=st.colors[1])
            ax0.plot([t, t], [cR, cT], lw=1, c=st.colors[1], zorder=1.)

            ax1.plot(t_isi, adap_isi, lw=1, c=st.colors[1])
            ax1.plot([t, t], [adap, adap_after], c=st.colors[1], zorder=1.0)

            ax2.plot(t_isi, jpuff_isi, lw=1, c=st.colors[1])
            t_isi.clear()
            ca_isi.clear()
            jpuff_isi.clear()
            adap_isi.clear()

    #for i in range(nrISIs - 1):
    #    x_left = spike_times[i]
    #    x_right = spike_times[i + 1]
    #    dx = spike_times[i + 1] - spike_times[i]
    #    ax0.text(x_left + dx / 2, 1.15, f"$T_{i + 1}$", ha="center", va="center", clip_on=False)

    #    ax0.arrow(x_left + 0.05 * dx, 1.05, 0.9 * dx, 0, fc="k", length_includes_head=True, head_width=0.05,
    #              head_length=5.0, lw=0.5,
    #              clip_on=False)
    #    ax0.arrow(x_right - 0.05 * dx, 1.05, -0.9 * dx, 0, fc="k", length_includes_head=True, head_width=0.05,
    #              head_length=5.0, lw=0.5,
    #              clip_on=False)

    # ax0.set_xlabel("$t$ / s")
    ax0.set_ylabel(r"$c_i$")
    ax0.set_ylim([0.8 * cR, 1.2 * cT])
    ax0.set_yticks([cR, cT])
    ax0.set_yticklabels(["$c_R$", "$c_T$"])
    ax0.set_xticklabels([])
    ax0.set_xlim([0, spike_times[nrISIs - 1] + 20])

    ax1.set_ylabel(r"$\hat{c}_{er}$")
    ax1.set_xticklabels([])
    ax1.set_xlim([0, spike_times[nrISIs - 1] + 20])

    ax2.set_xlabel("$t$ / s")
    ax2.set_xlim([0, spike_times[nrISIs - 1] + 20])
    ax2.set_ylabel(r"$j_{\rm puff}$")

    nr_ISIs = len(ISIs)
    index_ISIs = np.arange(nr_ISIs)
    popt, pcov = curve_fit(transient_func, index_ISIs, ISIs, p0=(100, 150, 2))
    print(popt)
    ISI_fit = popt[0] * np.exp(-index_ISIs / popt[2]) + popt[1] * (1. - np.exp(-index_ISIs / popt[2]))

    ax3.set_xlim([0, 50])
    ax3.set_ylim([0, 100])
    ax3.set_xlabel("$i$")
    ax3.set_ylabel("$T_i$ /s")
    cas = [ca for ca in cas if ca != 1]
    ax3.scatter(range(len(ISIs))[:50], ISIs[:50], fc="w", ec=st.colors[1], s=20, zorder=3)
    ax3.plot(index_ISIs, ISI_fit, lw=1, c="k")
    ax3.axhline(popt[0], ls=":", lw=1, c="k")
    ax3.axhline(popt[1], ls=":", lw=1, c="k")
    ax3.text(25, popt[0] * 1.2, "$T_0$", ha="center")
    ax3.text(25, popt[1] * 1.2, "$T_\infty$", ha="center")

    ax4.set_xlabel("$T_i$ / s")
    ax4.set_ylabel("$P_0(T_i)$")
    ax4.hist(ISIs, bins=20, color=st.colors[1], density=True, alpha=0.6)
    ts_inv_gau = np.linspace(0, 1.75 * mean_isi, 1001)
    inv_gaus = []
    for t in ts_inv_gau:
        p = np.sqrt(mean_isi / (2 * np.pi * np.power(cv_isi, 2) * (t ** 3))) * np.exp(
            -(t - mean_isi) ** 2 / (2 * mean_isi * np.power(cv_isi, 2) * t))
        inv_gaus.append(p)
    ax4.plot(ts_inv_gau, inv_gaus, ls="--", color="k", label="Inv.\ Gaussian")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig1.pdf", transparent=True)
    plt.show()
