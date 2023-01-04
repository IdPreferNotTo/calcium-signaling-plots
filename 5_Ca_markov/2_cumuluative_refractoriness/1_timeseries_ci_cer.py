import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(8, 4))
    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    st.remove_top_right_axis([ax1, ax2])
    set = 2

    if set == 1:
        tau = 1.0
        p = 0.060
        tau_er = 100
        eps_er = 0.1
    if set == 2:
        tau = 5.0
        p = 0.015
        tau_er = 100
        eps_er = 0.1

    home = os.path.expanduser("~")
    print("Load Data...")
    data_ci = df.load_traces_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
    data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
    mean_isi = np.mean(data_isi)
    ts, cas, jpuffs, adaps = np.transpose(data_ci)
    print("done")

    cer_top = []
    cer_bot = []
    for ca, jp, ap in zip(cas, jpuffs, adaps):
        if ca == 0.5 and jp == 0.:
            cer_top.append(ap)
            cer_bot.append((1-eps_er)*ap)

    mean_cer = np.mean(adaps[10_000:])
    mean_cer_top = np.mean(cer_top[10:])
    mean_cer_bot = np.mean(cer_bot[10:])
    sample_size = len(cer_top[10:])
    std = np.std(cer_top[10:])

    mean_cer_top_theory = 1 - eps_er*np.exp(-mean_isi/tau_er)/(1 - (1-eps_er)*np.exp(-mean_isi/tau_er))
    mean_cer_bot_theory = 1 - eps_er/(1 - (1-eps_er)*np.exp(-mean_isi/tau_er))
    mean_cer_theory = 1/(1 + (tau_er/mean_isi)*(eps_er/(1 - eps_er/2)))

    print(1 - mean_cer_top_theory, (1 - mean_cer_bot_theory)*np.exp(-mean_isi/tau_er))
    print(mean_cer_top_theory*(1.-eps_er), mean_cer_bot_theory)

    spike_times = []
    cR = 0.2
    cT = 0.5
    ts_plot = []
    cas_plot = []
    as_plot = []
    jp_plot = []
    n_spikes = 0
    for t, c, jp, a in zip(ts, cas, jpuffs, adaps):
        if t >= 2_000:
            break
        ts_plot.append(t)
        as_plot.append(a)
        cas_plot.append(c)
        jp_plot.append(jp)
        if c == 0.5 and jp == 0:
            spike_times.append(t)
            ts_plot.append(t)
            cas_plot.append(0.5 + a / 2)
            as_plot.append(a)
            jp_plot.append(jp)

            ts_plot.append(t)
            cas_plot.append(0.2)
            as_plot.append(a * (1 - eps_er))
            jp_plot.append(jp)
            n_spikes += 1

    ax1.set_ylabel(r"$c_{er}$")
    ax1.set_xlim([1000, 2000])
    ax1.plot(ts_plot, as_plot)
    ax1.axhline(mean_cer_bot, c="C7", ls="--", lw=2)
    ax1.axhline(mean_cer_top, c="C7", ls="--", lw=2)
    ax1.axhline(mean_cer, c="C7", ls="--", lw=2)

    ax1.axhline(mean_cer_bot_theory, c="C3")
    ax1.axhline(mean_cer_top_theory, c="C3")
    ax1.axhline(mean_cer_theory, c="C1")

    ax2.set_ylabel(r"$c_i$")
    ax2.set_ylim([0.8*cR, 2.1*cT])
    ax2.set_yticks([cR, cT])
    ax2.set_yticklabels(["$c_R$", "$c_T$"])
    ax2.set_xlim([1000, 2000])
    ax2.plot(ts_plot, cas_plot)
    plt.show()
