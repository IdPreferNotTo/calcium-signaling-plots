import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os
from scipy.optimize import curve_fit

import styles as st
import functions as fc

def func(x, a):
    return a/(np.sqrt(1 + 0.175*x))

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    fig = plt.figure(tight_layout=True, figsize=(6, 2.5))
    gs = gridspec.GridSpec(nrows=1, ncols=2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    st.remove_top_right_axis([ax1, ax2])
    colors = st.Colors()

    ax1.text(0.10, 0.95, r"A", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.10, 0.95, r"B", fontsize=13, transform=ax2.transAxes, va='top')

    home = os.path.expanduser("~")


    ax1.set_ylabel(r"$\langle T \rangle$")
    ax1.set_xlabel(r"$b_T$")
    #ax5.set_ylim([0, 100])

    ax2.set_ylabel(r"$CV_T$")
    ax2.set_xlabel(r"$b_T$")
    ax2.set_ylim([0, 0.5])

    tau = 5.62
    j = 0.0126

    bTs = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    #bTs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    mean_isis_buffer_markov = []
    cv_isis_buffer_markov = []
    folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_markov_buffer/out/kp10km50/"
    file_spikes = f"spike_times_markov_buffer_bt0.00_ip1.00_tau5.62e+00_j1.26e-02_K10_5.dat"
    data_isis = np.loadtxt(folder + file_spikes)
    mean0 = np.mean(data_isis)
    mean_bT = []
    for bt in bTs:
        beta = 1 + 5*bt/np.power(5 + 0.35, 2)
        mean_bT.append(mean0 * beta)

    for bT in bTs:
        folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_markov_buffer/out/kp10km50/"
        file_spikes = f"spike_times_markov_buffer_bt{bT:.2f}_ip1.00_tau5.62e+00_j1.26e-02_K10_5.dat"
        data_isis = np.loadtxt(folder + file_spikes)
        mean = np.mean(data_isis)
        cv = np.std(data_isis) / mean
        mean_isis_buffer_markov.append(mean)
        cv_isis_buffer_markov.append(cv)

    ax1.plot(bTs, mean_bT, ls="--", color=colors.palette[5], label="Theory")
    ax1.plot(bTs, mean_isis_buffer_markov, color=st.colors[0], zorder=1)
    ax2.plot(bTs, cv_isis_buffer_markov,color=st.colors[0], zorder=1, label="Two-component")
    cb0s = []
    p_cb0s = []
    mean_isis_buffer = []
    cv_isis_buffer = []
    for bT in bTs:
        folder_langevin = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin_buffer/out/kp10km50/"
        file_spikes_langevin = f"spike_times_langevin_buffer_bt{bT:.2f}_ip1.00_tau5.62e+00_j1.26e-02_K10_5.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_langevin = np.mean(data_isis_langevin)
        cv_langevin = np.std(data_isis_langevin) / mean_langevin
        mean_isis_buffer.append(mean_langevin)
        cv_isis_buffer.append(cv_langevin)

    popt, pcov = curve_fit(func, bTs, cv_isis_buffer)
    print(popt)
    c1 = popt[0]

    ax1.plot()
    ax1.scatter(bTs, mean_isis_buffer, s=20, fc="w", ec=st.colors[0], zorder=2, label="Sim.")
    ax2.scatter(bTs, cv_isis_buffer, s=20, fc="w", ec=st.colors[0], zorder=2, label="Langevin")
    cv_theory = []
    for bT in bTs:
        cv_theory.append(0.207992/np.sqrt(1 + 0.175*bT))
    ax2.plot(bTs, cv_theory, ls="--", color=colors.palette[5], label="Theory", zorder=1)
    ax2.legend(fancybox=False, fontsize=9, loc=1)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig11.pdf", transparent=True)
    plt.show()