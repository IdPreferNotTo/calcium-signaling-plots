import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os
from functions import *


if __name__ == "__main__":
    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, (4/3) * 6 / 2))
    gs = gridspec.GridSpec(4, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1:3])
    ax3 = fig.add_subplot(gs[3])
    remove_top_right_axis([ax1, ax2, ax3])

    tau = 2.81
    j = 0.0728
    taua = 100
    ampa = 0.2
    home = os.path.expanduser("~")

    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike =  f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    isis  = np.loadtxt(home + folder + file_spike)
    mean_isi = np.mean(isis)
    cv_isi = np.std(isis)/mean_isi
    print(mean_isi, cv_isi)

    ts, cas, jpuffs, adaps = np.transpose(data)

    spike_times = []
    cR = 0.33
    cT = 1.0
    t_isi = []
    ca_isi = []
    jpuff_isi = []
    adap_isi= []
    for t, ca, jpuff, adap, adap_after in zip(ts[:-1], cas[:-1], jpuffs[:-1], adaps[:-1], adaps[1:]):
        t_isi.append(t)
        ca_isi.append(ca)
        jpuff_isi.append(jpuff)
        adap_isi.append(adap)
        if ca == 1:
            spike_times.append(t)
            ax1.plot(t_isi, adap_isi, c="C0")
            ax1.plot([t, t], [adap, adap_after], ls=":", c="C0")

            ax2.arrow(x=t, y=1, dx=0, dy=1, color="C0", length_includes_head = True, head_width = 7.5, head_length=0.1)
            ax2.plot(t_isi, ca_isi, c="C0")
            ax2.plot([t, t], [cR, cT], ls=":", c="C0")

            ax3.plot(t_isi, jpuff_isi, c="C0")
            t_isi.clear()
            ca_isi.clear()
            jpuff_isi.clear()
            adap_isi.clear()


    for i in range(7):
        x_left = spike_times[i]
        x_right = spike_times[i+1]
        dx = spike_times[i+1] - spike_times[i]
        ax2.text(x_left + dx/2, 1.6, f"$T_{i+1}$", ha="center", va="center", clip_on=False)

        ax2.arrow(x_left + 0.05*dx, 1.5, 0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=15.0, lw=0.5,
                clip_on=False)
        ax2.arrow(x_right -0.05*dx, 1.5, -0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=15.0, lw=0.5,
                clip_on=False)

    ax1.set_ylabel(r"[Ca\textsuperscript{2+}]$_{\rm ER}$ [a.u.]")
    ax1.set_xlim([0, 500])
    ax1.set_ylim([0.5, 1])

    ax2.set_ylabel(r"[Ca\textsuperscript{2+}]$_{i}$ [a.u.]")
    ax2.set_ylim([0.8*cR, 2.1*cT])
    ax2.set_yticks([cR, cT])
    ax2.set_yticklabels(["$c_R$", "$c_T$"])
    ax2.set_xlim([0, 500])

    ax3.set_xlabel("$t$ [s]")
    ax3.set_xlim([0, 500])
    ax3.set_ylabel(r"$J_{\rm puff}$ [a.u.]")
    plt.savefig(home + f"/Data/Calcium/Plots/7_markov_ca_adap_timeseries_tau{tau:.2e}j{j:.2e}.pdf", transparent=True)

    plt.show()
