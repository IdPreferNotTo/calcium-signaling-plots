import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(3, 1)
    ax = fig.add_subplot(gs[0:2])
    ax2 = fig.add_subplot(gs[2])
    st.remove_top_right_axis([ax, ax2])

    # Premium parameters should definitely be remembered
    tau = 5
    j = 0.015
    data = df.load_traces_markov(tau, j, cer=False)
    isis  = df.load_spike_times_markov(tau, j, cer=False)
    mean_isi = np.mean(isis)
    cv_isi = np.std(isis)/mean_isi
    print(mean_isi, cv_isi)

    ts, cas, jpuffs, adaps = np.transpose(data)

    #ax.plot(ts, cas)

    spike_times = []
    cR = 0.20
    cT = 0.50
    t_plot = []
    ca_plot = []
    jp_plot = []
    for t, ca, jp in zip(ts, cas, jpuffs):
        if t > 250:
            break
        t_plot.append(t)
        ca_plot.append(ca)
        jp_plot.append(jp)
        if ca == cT and jp == 0:
            spike_times.append(t)

            t_plot.append(t)
            ca_plot.append(ca + 1)
            jp_plot.append(jp)
            t_plot.append(t)
            ca_plot.append(cR)
            jp_plot.append(jp)

    ax.plot(t_plot, ca_plot, lw=1, color=st.colors[0], zorder=1)
    ax2.plot(t_plot, jp_plot, c=st.colors[0])

    for i in range(4):
        if i == 0:
            x_left = 0
        else:
            x_left = spike_times[i - 1]
        x_right = spike_times[i]
        dx = x_right - x_left
        ax.text(x_left + dx/2, 0.9, f"$T_{i+1}$", ha="center", va="center", clip_on=False)
        ax.annotate("" , xy=(x_right, 0.8), xytext=(x_left, 0.8),  arrowprops=dict(arrowstyle='<->', connectionstyle="arc3", color="k", lw=1))

    ax.set_xlim([0, 250])
    ax2.set_xlabel("$t$")
    ax2.set_xlim([0, 250])
    ax2.set_ylabel(r"$j_{\rm puff}$")
    ax2.set_ylim([-0.05, 0.75])
    ax.set_ylabel(r"$c_{\rm i}$")
    ax.set_ylim([0.8*cR, 2.0*cT])
    ax.set_yticks([cR, cT])
    ax.set_yticklabels(["$c_R$", "$c_T$"])
    plt.show()
