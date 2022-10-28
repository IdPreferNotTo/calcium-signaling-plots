import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2)
    ax1a = fig.add_subplot(gs[0, 0])
    ax1b = fig.add_subplot(gs[0, 1])
    ax2a = fig.add_subplot(gs[1, 0])
    ax2b = fig.add_subplot(gs[1, 1])
    axis = [ax1a, ax1b, ax2a, ax2b]
    st.remove_top_right_axis(axis)
    ax1a.set_ylabel(r"$p_{ISI}(t)$")
    ax1a.set_xlabel(r"$t$")
    ax2a.set_ylabel(r"$p_{ISI}(t)$")
    ax2a.set_xlabel(r"$t$")
    # Parameters
    tau = 5
    p = 0.015
    tau_er = 501
    eps_ers = [0.0398, 0.398]

    pbif = fc.pbif(tau, K=10, N=5, M=3, s=1.0)
    cmap_cividis = plt.get_cmap("YlGnBu", 10)
    for eps_er, axa, axb in zip(eps_ers, [ax1a, ax2a], [ax1b, ax2b]):
        cer_crit = pbif/p
        data_ci = df.load_traces_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
        data_exit = [data for data in data_ci if data[1] == 0.5 and data[2] == 0]

        data_exit_excitable = [data for data in data_ci if data[1] == 0.5 and data[2] == 0 and data[3] < cer_crit]
        data_exit_meandriven = [data for data in data_ci if data[1] == 0.5 and data[2] == 0 and data[3] > cer_crit]
        data_exit = np.asarray(data_exit)


        N, bins, patches = axa.hist(data_exit[5:,3], density=True, alpha=.75, color=st.colors[0], bins=50)
        for idx, bin in enumerate(bins):
            if bin > cer_crit:
                break
        for i in range(idx):
            patches[i].set_facecolor(st.colors[5])
        for i in range(idx, len(patches)):
            patches[i].set_facecolor(st.colors[1])
        data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)

        times, cis, jps, cers = np.transpose(data_exit)
        isi_excitable = []
        isi_mean_driven = []
        for isi, cer in zip(data_isi, cers):
            if cer > cer_crit:
                isi_mean_driven.append(isi)
            else:
                isi_excitable.append(isi)

        mean = np.mean(data_isi)
        std = np.std(data_isi)
        cv = std/mean
        cv2 = cv**2
        ts = np.linspace(mean-3*std, mean+3*std, 501)
        axb.hist([isi_mean_driven, isi_excitable], bins=50, color=[st.colors[1], st.colors[5]], stacked=True, alpha=0.75, density=True)

    plt.show()