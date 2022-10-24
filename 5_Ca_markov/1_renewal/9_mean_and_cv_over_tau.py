import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 2.4))
    gs = gridspec.GridSpec(nrows=1, ncols=2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    axis = [ax1, ax2]
    st.remove_top_right_axis(axis)

    ax1.text(0.1, 1.00, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 1.00, r"B", fontsize=11, transform=ax2.transAxes, va='top')

    K = 10
    N = 5
    M = 3
    tau = 1

    MEANs = []
    CVs = []
    js = np.logspace(-2, 0, 101)
    for j in js:
        N += 1
        file_str = home + f"/Data/calcium_spikes_markov/Data_no_adap/spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K{K:d}_5.dat"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ISIs = np.loadtxt(file_str, unpack=True)
            if ISIs.size == 0:
                T = np.infty
                CV = 1.
                n = ISIs.size
            else:
                T = np.mean(ISIs)
                CV = np.std(ISIs) / np.mean(ISIs)
                n = ISIs.size
        MEANs.append(T)
        CVs.append(CV)
    ax1.set_yscale("log")
    ax1.set_ylim([0.1, 100])
    ax1.set_xlim([0, 1])
    ax1.plot(js, MEANs)

    ax2.plot(js, CVs)
    ax2.set_xlim([0, 1])
    plt.show()