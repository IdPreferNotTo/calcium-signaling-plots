import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os
from scipy.stats import moment
from scipy.optimize import curve_fit

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
        err_means = []
        cvs = []
        err_cvs = []
        folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_markov_buffer/out/"
        file_spikes = f"spike_times_markov_buffer_bt0.00_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        mean0 = np.mean(data_isis)
        cv0 = np.std(data_isis)/mean0
        print(cv0)
        mean_bT = []
        for bt in bTs:
            beta = 1 + 5*bt/np.power(5 + 0.35, 2)
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

            cvs_tmp = []
            means_tmp = []
            n = 100
            for isi_sample in [data_isis[i:i + n] for i in range(0, len(data_isis), n)]:
                mean_tmp = np.mean(isi_sample)
                cv_tmp = np.std(isi_sample) / mean_tmp
                means_tmp.append(mean_tmp)
                cvs_tmp.append(cv_tmp)
            #err_mean = np.sqrt(var/N)
            #err_cv = np.sqrt((moment(data_isis, moment=4) - var**2)/N)
            #err_means.append(err_mean)
            #err_cvs.append(err_cv)
            err_means.append(np.std(means_tmp))
            err_cvs.append(np.std(cvs_tmp))
        ax_mean.plot(bTs, mean_bT, color=colors.palette[5], label="Theory")
        ax_mean.errorbar(bTs, means, yerr=err_means, linestyle="none", marker="o", ms=4.5, mfc="w", mec=st.colors[0], zorder=2, label="Two-component \n + fast buffer")
        ax_cv.errorbar(bTs, cvs, yerr=err_cvs, linestyle="none", marker="o", ms=4.5, mfc="w", mec=st.colors[0], zorder=2, label="Two-component \n + fast buffer")
        if tau == 5.0:
            ax_mean.set_ylim([0, 250])
            ax_cv.set_ylim([0, 0.35])
            cv_theory = []
            for bT in np.linspace(0, 100, 100):
                cv_theory.append(0.176768/np.sqrt(1 + 0.175*bT))
            ax_cv.plot(np.linspace(0, 100, 100), cv_theory, color=colors.palette[5], label="Theory", zorder=1)
        if tau == 1.0:
            ax_mean.set_ylim([0, 700])
            ax_cv.set_ylim([0, 0.7])
            cv_theory = []
            for bT in np.linspace(0, 100, 100):
                cv_theory.append(cv0 / np.sqrt(1 + 0.175 * bT))
            ax_cv.plot(np.linspace(0, 100, 100), cv_theory, color=colors.palette[5], label="Theory", zorder=1)
    ax2.legend(fancybox=False, fontsize=9, loc=1)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig11.pdf", transparent=True)
    plt.show()