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
    j = 0.025
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_markov_buffer/out/"
    file_ca = f"ca_markov_bt50.00_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
    data = np.loadtxt(folder + file_ca)

    ts, cas, cbs, jpuffs, adaps = np.transpose(data)

    #ax.plot(ts, cas)

    spike_times = []
    cR = 0.20
    cT = 0.50
    t_plot = []
    ca_plot = []
    jp_plot = []
    for t, ca, cb, jp in zip(ts, cas, cbs, jpuffs):
        if t > 1000:
            break
        t_plot.append(t)
        ca_plot.append(ca)
        jp_plot.append(cb)
        if ca == cT and jp == 0:
            spike_times.append(t)

            t_plot.append(t)
            ca_plot.append(ca + 1)
            jp_plot.append(cb)
            t_plot.append(t)
            ca_plot.append(cR)
            jp_plot.append(cb)

    ax.plot(t_plot, ca_plot, lw=1, color=st.colors[0], zorder=1)
    ax2.plot(t_plot, jp_plot, lw=1, c=st.colors[0])

    for i in range(4):
        if i == 0:
            x_left = 0
        else:
            x_left = spike_times[i - 1]
        x_right = spike_times[i]
        dx = x_right - x_left
        ax.text(x_left + dx/2, 0.9, f"$T_{i+1}$", ha="center", va="center", clip_on=False)
        ax.annotate("" , xy=(x_right, 0.8), xytext=(x_left, 0.8),  arrowprops=dict(arrowstyle='<->', connectionstyle="arc3", color="k", lw=1))

    ax.set_xlim([0, 1000])
    ax2.set_xlabel("$t$")
    ax2.set_xlim([0, 1000])
    ax2.set_ylabel(r"$c_b$")
    ax2.set_ylim([0.0, 5.00])
    ax.set_ylabel(r"$c_{\rm i}$")
    ax.set_ylim([0.8*cR, 2.0*cT])
    ax.set_yticks([cR, cT])
    ax.set_yticklabels(["$c_R$", "$c_T$"])
    plt.savefig(home + f"/Plots/timeseries_ci_cb.pdf")
    plt.show()
