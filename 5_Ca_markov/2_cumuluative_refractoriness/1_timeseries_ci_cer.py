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
    gs = gridspec.GridSpec(4, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1:3, 0])
    ax3 = fig.add_subplot(gs[3, 0])
    ax4 = fig.add_subplot(gs[0:2, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    tau = 5.0
    p = 0.015
    tau_er = 100
    eps_er = 0.1
    home = os.path.expanduser("~")
    print("Load Data...")
    data_ci = df.load_traces_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
    data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)

    ts, cas, jpuffs, adaps = np.transpose(data_ci)
    print("done")

    spike_times = []
    cR = 0.2
    cT = 0.5
    ts_plot = []
    cas_plot = []
    as_plot = []
    jp_plot = []
    n_spikes = 0
    for t, c, jp, a in zip(ts, cas, jpuffs, adaps):
        if t >= 1_000:
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
    ax1.set_xlim([0, 500])
    ax1.set_ylim([0.5, 1])
    ax1.plot(ts_plot, as_plot)

    ax2.set_ylabel(r"$c_i$")
    ax2.set_ylim([0.8*cR, 2.1*cT])
    ax2.set_yticks([cR, cT])
    ax2.set_yticklabels(["$c_R$", "$c_T$"])
    ax2.set_xlim([0, 500])
    ax2.plot(ts_plot, cas_plot)

    ax3.set_xlabel("$t$ / s")
    ax3.set_xlim([0, 500])
    ax3.set_ylabel(r"$j_{\rm puff}$")
    ax3.plot(ts_plot, jp_plot)

    ax4.hist(data_isi, bins=50, density=True)
    mean_isi = np.mean(data_isi)
    cv_isi = np.std(data_isi)/mean_isi

    ts = np.linspace(0, 2*mean_isi, 100)
    density_inverse_gaus = fc.inverse_gaussian_dist(ts, mean_isi, cv_isi**2)
    ax4.plot(ts, density_inverse_gaus)

    plt.show()
