import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 3))
    gs = gridspec.GridSpec(2, 2)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3])

    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_no_adap/"
    file = "ca_markov_ip1.00_tau1.05e+01_j1.10e-02_N10_0.dat"
    file_spikes = "spike_times_markov_ip1.00_tau1.05e+01_j1.10e-02_N10_0.dat"
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
    adap_isi= []
    nrISIs = 0
    for t, ca, jpuff, adap, adap_after in zip(ts[:-1], cas[:-1], jpuffs[:-1], adaps[:-1], adaps[1:]):
        t_isi.append(t)
        ca_isi.append(ca)
        jpuff_isi.append(jpuff)
        adap_isi.append(adap)
        if ca == 1 and nrISIs <= 4:
            nrISIs += 1
            spike_times.append(t)
            ax0.plot([t, t], [adap, adap_after], c=st.colors[0], zorder=1.0)

            ax0.plot(t_isi, ca_isi, lw=1, c=st.colors[1])
            ax0.plot([t, t], [cR, cT], lw=1, c=st.colors[1], zorder=1.)

            ax2.plot(t_isi, jpuff_isi, lw=1, c=st.colors[1])
            t_isi.clear()
            ca_isi.clear()
            jpuff_isi.clear()
            adap_isi.clear()


    for i in range(nrISIs-1):
        x_left = spike_times[i]
        x_right = spike_times[i+1]
        dx = spike_times[i+1] - spike_times[i]
        ax0.text(x_left + dx/2, 1.15, f"$T_{i+1}$", ha="center", va="center", clip_on=False)

        ax0.arrow(x_left + 0.05*dx, 1.05, 0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                clip_on=False)
        ax0.arrow(x_right -0.05*dx, 1.05, -0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                clip_on=False)

    #ax0.set_xlabel("$t$ / s")
    ax0.set_ylabel(r"$c_i$")
    ax0.set_ylim([0.8*cR, 1.5*cT])
    ax0.set_yticks([cR, cT])
    ax0.set_yticklabels(["$c_R$", "$c_T$"])
    ax0.set_xticklabels([])
    ax0.set_xlim([0, spike_times[nrISIs-1] + 20])

    ax2.set_xlabel("$t$ / s")
    ax2.set_xlim([0, spike_times[nrISIs-1] + 20])
    ax2.set_ylabel(r"$j_{\rm puff}$")

    file_p0_theory = home + f"/Data/Calcium/data/p0_ip1.00_tau1.05e+01_j1.10e-02_N10.dat"
    data = np.loadtxt(file_p0_theory)
    ca_theory, p0_theory = np.transpose(data)
    ax1.set_xlabel("$c_i$")
    ax1.set_ylabel("$P_0(c_i)$")
    cas = [ca for ca in cas if ca != 1]
    ax1.hist(cas, bins=25, density=True, alpha=.6, color=st.colors[1])
    ax1.set_xticks([0.33, 1.0])
    ax1.set_xticklabels(["$c_R$", "$c_T$"])
    ax1.plot(ca_theory[::100], p0_theory[::100], color=st.colors[0])

    ax3.set_xlabel("$T$ / s")
    ax3.set_ylabel("$P_0(T)$")
    ax3.hist(ISIs, bins=20, color=st.colors[1], density=True, alpha=0.6)
    ts_inv_gau = np.linspace(0, 1.75*mean_isi, 1001)
    inv_gaus = []
    for t in ts_inv_gau:
        p = np.sqrt(mean_isi / (2 * np.pi * np.power(cv_isi, 2) * (t ** 3))) * np.exp(
            -(t - mean_isi) ** 2 / (2 * mean_isi * np.power(cv_isi, 2) * t))
        inv_gaus.append(p)
    ax3.plot(ts_inv_gau, inv_gaus, ls="--", color="k", label="Inv.\ Gaussian")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig5.png",transparent=True)
    plt.show()
