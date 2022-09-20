import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as stats
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    # Set up plot style
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3, 4))
    gs = gridspec.GridSpec(nrows=1, ncols=1)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    st.remove_top_right_axis([ax1])

    ax1.set_ylabel("$\Delta c_{er}$")
    ax1.set_xlabel("$t$ / s")
    ax1.set_ylim([0.5, 1])

    ax2.set_ylabel("$T_i$")
    ax2.set_xlabel("$i$")
    ax2.set_ylim([0.5, 1])
    ca_r = 0.33
    ca_t = 1.
    n = 5
    m = 4
    N = 10
    ip3 = 1.0

    #Plot regular spike train
    tau_i = 10.
    p = 0.015
    eps = 0.1
    tau_a = 500
    home = os.path.expanduser("~")
    folder = home + "/Data/calcium_spikes_markov/Data_adap/"
    file = f"ca_markov_ip1.00_taua{tau_a:.2e}_ampa{eps:.2e}_tau{tau_i:.2e}_j{p:.2e}_N10_0.dat"
    file_isi = f"/spike_times_markov_ip1.00_taua{tau_a:.2e}_ampa{eps:.2e}_tau{tau_i:.2e}_j{p:.2e}_N10_0.dat"
    data = np.loadtxt(folder + file)
    isis = np.loadtxt(folder + file_isi)
    t_s, cai_s, jpuff_s, cer_s = np.transpose(data)
    T_infty = np.mean(isis[100:])

    a_infty_min = eps/(1 - (1-eps)*np.exp(-T_infty/tau_a))
    a_infty_max = eps*np.exp(-T_infty/tau_a)/(1 - (1-eps)*np.exp(-T_infty/tau_a))
    ax1.axhline(1 - a_infty_min, lw=1, ls=":", color="C7", zorder=2)
    ax1.axhline((1 - a_infty_max), lw=1, ls=":", color="C7", zorder=2)

    spike_times = []
    cer_i_s = []
    count = 0
    max_spikes = 7
    max_i = 0
    while count < max_spikes:
        if cai_s[max_i] == 1:
            count += 1

            cer_i_s.append(cer_s[max_i])
            spike_times.append(t_s[max_i])
            ax1.text(t_s[max_i], cer_s[max_i] * (1 - eps) - 0.025, f"$\Delta c_{count:d}$", fontsize=11, va='top', ha="center", zorder=1)
            ax1.scatter(t_s[max_i], cer_s[max_i] * (1 - eps), ec=st.colors[0], fc="w", zorder=2)
            ax1.scatter(t_s[max_i], cer_s[max_i], marker="s", s=15, ec=st.colors[0], fc="w", zorder=2)
        max_i += 1


    ax1.plot(t_s[0:max_i], cer_s[0:max_i], color=st.colors[0], zorder=1)
    ax1.set_ylim([0.5, 1.1])
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig9.png")
    plt.show()