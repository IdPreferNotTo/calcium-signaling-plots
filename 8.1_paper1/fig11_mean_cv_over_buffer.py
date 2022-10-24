import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os

import styles as st
import functions as fc

def func(x, a):
    return a/(np.sqrt(1 + 0.175*x))

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    fig = plt.figure(tight_layout=True, figsize=(6, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])
    colors = st.Colors()

    ax1.text(0.10, 0.95, r"A$_1$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.10, 0.95, r"A$_2$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.10, 0.95, r"B$_1$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.10, 0.95, r"B$_2$", fontsize=13, transform=ax4.transAxes, va='top')

    home = os.path.expanduser("~")

    axis_mean = [ax1, ax3]
    axis_cv = [ax2, ax4]


    bTs = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    taus = [5, 1]
    js = [0.025, 0.064]
    for ax_mean, ax_cv, tau, j in zip(axis_mean, axis_cv, taus, js):
        ax_mean.set_ylabel(r"$\langle T \rangle$")
        ax_mean.set_xlabel(r"$b_T$")
        # ax5.set_ylim([0, 100])

        ax_cv.set_ylabel(r"$CV_T$")
        ax_cv.set_xlabel(r"$b_T$")
        means = []
        means_langevin = []
        err_means = []
        cvs = []
        cvs_langevin = []
        err_cvs = []
        folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_markov_buffer/out/"
        file_spikes = f"spike_times_markov_buffer_bt0.00_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        mean0 = np.mean(data_isis)
        cv0 = np.std(data_isis)/mean0
        print(cv0)
        mean_bT = []
        for bt in bTs:
            beta = 1 + 5*bt/np.power(5, 2)
            mean_bT.append(mean0 * beta)

        for bT in bTs:
            folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_markov_buffer/out/"
            file_spikes = f"spike_times_markov_buffer_bt{bT:.2f}_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
            data_isis = np.loadtxt(folder + file_spikes)
            N = np.size(data_isis)
            mean = np.mean(data_isis)
            var = np.var(data_isis)
            cv = np.std(data_isis) / mean
            means.append(mean)
            cvs.append(cv)

            folder_langevin = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin_buffer/out/"
            file_spikes_langevin = f"spike_times_langevin_buffer_bt{bT:.2f}_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_0.dat"
            data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
            N_langevin = np.size(data_isis_langevin)
            mean_langevin = np.mean(data_isis_langevin)
            var_langevin = np.var(data_isis_langevin)
            cv_langevin = np.std(data_isis_langevin) / mean_langevin
            means_langevin.append(mean_langevin)
            cvs_langevin.append(cv_langevin)

        ax_mean.scatter(bTs, means, s=20, fc="w", ec=st.colors[0], zorder=2)
        ax_mean.plot(bTs, means_langevin, color=colors.palette[0], zorder=1, label="Langevin + buffer")
        ax_mean.plot(bTs, mean_bT, color=colors.palette[5], zorder=0, label="Theory")
        ax_cv.scatter(bTs, cvs, s=20, fc="w", ec=st.colors[0], zorder=2, label="Two-component \n + buffer")
        ax_cv.plot(bTs, cvs_langevin, color=colors.palette[0], zorder=1, label="Langevin \n + buffer")
        if tau == 5.0:
            ax_mean.set_ylim([0, 250])
            ax_cv.set_ylim([0, 0.35])
            cv_theory = []
            for bT in np.linspace(0, 100, 100):
                cv_theory.append(0.176768/np.sqrt(1 + 5*bT/np.power(5, 2)))
            ax_cv.plot(np.linspace(0, 100, 100), cv_theory, color=colors.palette[5], label="Theory", zorder=0)
        if tau == 1.0:
            ax_cv.set_ylim([0, 0.8])

    handles, labels = ax2.get_legend_handles_labels()
    order = [0, 1, 2]
    ax2.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fancybox=False, fontsize=9, loc=1)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig11.pdf", transparent=True)
    plt.show()