import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
from functions import *


if __name__ == "__main__":
    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(3, 1)
    ax = fig.add_subplot(gs[0:2])
    ax2 = fig.add_subplot(gs[2])
    remove_top_right_axis([ax, ax2])

    # Premium parameters should definitely be remembered
    tau = 2.81
    j = 0.0574
    home = os.path.expanduser("~")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/"
    file = f"ca_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike =  f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    isis  = np.loadtxt(home + folder + file_spike)
    mean_isi = np.mean(isis)
    cv_isi = np.std(isis)/mean_isi
    print(mean_isi, cv_isi)

    ts, cas, jpuffs, adaps = np.transpose(data)

    #ax.plot(ts, cas)

    spike_times = []
    cR = 0.33
    cT = 1.0
    t_isi = []
    ca_isi = []
    jpuff_isi = []
    for t, ca, jpuff in zip(ts, cas, jpuffs):
        t_isi.append(t)
        ca_isi.append(ca)
        jpuff_isi.append(jpuff)
        if ca == 1:
            spike_times.append(t)
            ax.arrow(x=t, y=1, dx=0, dy=1, color="C0", length_includes_head = True, head_width = 7.5, head_length=0.1)
            ax.plot(t_isi, ca_isi, c="C0")
            ax.plot([t, t], [cR, cT], ls=":", c="C0")

            ax2.plot(t_isi, jpuff_isi, c="C0")
            #ax2.plot([t, t], [cR, cT], ls=":", c="C0")
            t_isi.clear()
            ca_isi.clear()
            jpuff_isi.clear()


    for i in range(7):
        x_left = spike_times[i]
        x_right = spike_times[i + 1]
        dx = spike_times[i + 1] - spike_times[i]
        ax.text(x_left + dx/2, 1.6, f"$T_{i+1}$", ha="center", va="center", clip_on=False)

        ax.arrow(x_left + 0.05*dx, 1.5, 0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=15.0, lw=0.5,
                clip_on=False)
        ax.arrow(x_right -0.05*dx, 1.5, -0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=15.0, lw=0.5,
                clip_on=False)

    ax.set_xlim([0, 500])
    ax2.set_xlabel("$t$ [s]")
    ax2.set_xlim([0, 500])
    ax2.set_ylabel(r"$J_{\rm puff}$ [a.u.]")

    ax.set_ylabel(r"[Ca\textsuperscript{2+}] [a.u.]")
    ax.set_ylim([0.8*cR, 2.1*cT])
    ax.set_yticks([cR, cT])
    ax.set_yticklabels(["$c_R$", "$c_T$"])

    plt.show()
