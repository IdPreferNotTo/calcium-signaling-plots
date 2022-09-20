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
    fig = plt.figure(tight_layout=True, figsize=(9, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    gs_reg = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs[0, 0], hspace=0.1)
    gs_irr = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs[0, 1], hspace=0.1)
    gs_tra = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs[0, 2], hspace=0.1)

    ax_ca_reg = fig.add_subplot(gs_reg[0:2])
    ax_ca_irr = fig.add_subplot(gs_irr[0:2])
    ax_ca_tra = fig.add_subplot(gs_tra[0:2])
    ax_adap_reg = fig.add_subplot(gs_reg[2])
    ax_adap_irr = fig.add_subplot(gs_irr[2])
    ax_adap_tra = fig.add_subplot(gs_tra[2])
    ax_f_reg = fig.add_subplot(gs[1, 0])
    ax_f_irr = fig.add_subplot(gs[1, 1])
    ax_f_tra = fig.add_subplot(gs[1, 2])
    st.remove_top_right_axis([ax_ca_reg, ax_ca_irr, ax_ca_tra, ax_adap_reg, ax_adap_irr, ax_adap_tra, ax_f_reg, ax_f_irr, ax_f_tra])

    for ax in [ax_ca_reg, ax_ca_irr, ax_ca_tra]:
        ax.set_ylim([0.2, 4])
        ax.set_yticks([0.33, 1])
        ax.set_yticklabels(["$c_R$", "$c_T$"])
        ax.set_ylabel("$c_i$")
        ax.axhline(1, lw=1, ls=":", color="C7")

    for ax in [ax_adap_reg, ax_adap_irr, ax_adap_tra]:
        ax.set_ylabel("$\Delta c_{er}$")
        ax.set_xlabel("$t$ / s")
        ax.set_ylim([0.5, 1])

    for ax in [ax_f_reg, ax_f_irr, ax_f_tra]:
        ax.set_xlim([0.33, 1.00])
        ax.set_xlabel("$c_i$")
        ax.set_ylabel("$f(c_i ; \Delta c_{er})$")
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
    ax_adap_reg.axhline(1-a_infty_min, lw=1, ls=":",color="C7")
    ax_adap_reg.axhline((1-a_infty_max), lw=1, ls=":",color="C7")

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
        max_i += 1

    ax_ca_reg.plot(t_s[0:max_i], cai_s[0:max_i], color=st.colors[0])
    ax_adap_reg.plot(t_s[0:max_i], cer_s[0:max_i], color=st.colors[0])

    fs = []
    for cer_i in cer_i_s:
        xs = np.linspace(0.3, 1, 100)
        ca_drifts = []
        for x in xs:
            ca_drift = -(x - ca_r) / tau_i + cer_i * p * N * fc.mean_puff_single(x, n, m, ip3)
            ca_drifts.append(ca_drift)
        ax_f_reg.plot(xs, ca_drifts)
        fs.append(ca_drifts)


    #Plot irregular spike train
    tau_i = 0.2
    p = 0.5
    tau_a = 500
    eps = 0.1
    home = os.path.expanduser("~")
    print("Load Data...")
    folder = "/Data/calcium_spikes_markov/Data_adap/"
    file = f"ca_markov_ip1.00_taua{tau_a:.2e}_ampa{eps:.2e}_tau{tau_i:.2e}_j{p:.2e}_N10_0.dat"
    file_isi = f"/spike_times_markov_ip1.00_taua{tau_a:.2e}_ampa{eps:.2e}_tau{tau_i:.2e}_j{p:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    data_isis = np.loadtxt(home + folder + file_isi)
    t_s, cai_s, jpuff_s, cer_s = np.transpose(data)
    T_infty = np.mean(data_isis[100:])

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
        max_i += 1

    ax_ca_irr.plot(t_s[0:max_i], cai_s[0:max_i], color=st.colors[0])
    ax_adap_irr.plot(t_s[0:max_i], cer_s[0:max_i], color=st.colors[0])
    a_infty_min = eps/(1 - (1-eps)*np.exp(-T_infty/tau_a))
    a_infty_max = eps*np.exp(-T_infty/tau_a)/(1 - (1-eps)*np.exp(-T_infty/tau_a))
    ax_adap_irr.axhline(1-a_infty_min, lw=1, ls=":",color="C7")
    ax_adap_irr.axhline((1-a_infty_max), lw=1, ls=":",color="C7")

    fs = []
    for cer_i in cer_i_s:
        xs = np.linspace(0.3, 1, 100)
        ca_drifts = []
        for x in xs:
            ca_drift = -(x - ca_r) / tau_i + cer_i * p * N * fc.mean_puff_single(x, n, m, ip3)
            ca_drifts.append(ca_drift)
        ax_f_irr.plot(xs, ca_drifts)
        fs.append(ca_drifts)


    #Plot transient spike train
    #Plot irregular spike train
    tau_i = 10.
    p = 0.015
    tau_a = 10_000
    eps = 0.1
    home = os.path.expanduser("~")
    print("Load Data...")
    folder = "/Data/calcium_spikes_markov/Data_adap/"
    file = f"ca_markov_ip1.00_taua{tau_a:.2e}_ampa{eps:.2e}_tau{tau_i:.2e}_j{p:.2e}_N10_0.dat"
    file_isi = f"/spike_times_markov_ip1.00_taua{tau_a:.2e}_ampa{eps:.2e}_tau{tau_i:.2e}_j{p:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    data_isis = np.loadtxt(home + folder + file_isi)
    t_s, cai_s, jpuff_s, cer_s = np.transpose(data)
    T_infty = np.mean(data_isis[100:])

    spike_times = []
    cer_i_s = []
    count = 0
    max_spikes = 7
    max_i = 0
    for ci, cer in zip(cai_s, cer_s):
        if ci == 1:
            count += 1
            cer_i_s.append(cer)

    ax_ca_tra.plot(t_s, cai_s, color=st.colors[0])
    ax_adap_tra.plot(t_s, cer_s, color=st.colors[0])
    a_infty_min = eps/(1 - (1-eps)*np.exp(-T_infty/tau_a))
    a_infty_max = eps*np.exp(-T_infty/tau_a)/(1 - (1-eps)*np.exp(-T_infty/tau_a))
    ax_adap_tra.axhline(1-a_infty_min, lw=1, ls=":",color="C7")
    ax_adap_tra.axhline((1-a_infty_max), lw=1, ls=":",color="C7")

    fs = []
    for cer_i in cer_i_s:
        xs = np.linspace(0.3, 1, 100)
        ca_drifts = []
        for x in xs:
            ca_drift = -(x - ca_r) / tau_i + cer_i * p * N * fc.mean_puff_single(x, n, m, ip3)
            ca_drifts.append(ca_drift)
        ax_f_tra.plot(xs, ca_drifts)
        fs.append(ca_drifts)
    plt.show()