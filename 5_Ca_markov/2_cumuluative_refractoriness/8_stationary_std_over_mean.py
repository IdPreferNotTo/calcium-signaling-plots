import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axis_top = [ax1, ax2]
    axis_bot = [ax3, ax4]
    axis = axis_top + axis_bot
    st.remove_top_right_axis(axis)

    for ax in axis:
        ax.set_xlabel(r"$\langle T \rangle$")
        ax.set_ylabel(r"$\sigma_T$")
        #ax.set_xscale("log")
        #ax.set_yscale("log")
        ax.set_xlim([0, 250])
        ax.set_ylim([0, 100])

    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.1
    tau_er_fix = 100
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)

    for i, (tau, p) in enumerate(zip(taus, ps)):
        means = []
        stds = []
        p1s = []
        for tau_er in tau_ers:
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er_fix)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            stds.append(std)
            p1s.append(p1)

        if i ==0:
            ax1.scatter(means, stds, fc="w", ec=st.colors[i], s=20, zorder=3)
        else:
            ax3.scatter(means, stds, fc="w", ec=st.colors[i], s=20, zorder=3)

        means = []
        stds = []
        p1s = []
        for eps_er in eps_ers:
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er_fix, ampa=eps_er)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            stds.append(std)
            p1s.append(p1)

        if i == 0:
            ax2.scatter(means, stds, fc="w", ec=st.colors[i], s=20, zorder=3)
        else:
            ax4.scatter(means, stds, fc="w", ec=st.colors[i], s=20, zorder=3)

    plt.show()
