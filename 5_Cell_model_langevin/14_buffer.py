import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    fig = plt.figure(tight_layout=True, figsize=(6, 4.5))
    gs = gridspec.GridSpec(nrows=2, ncols=2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    home = os.path.expanduser("~")
    bTs = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    bTs = [10, 50]
    cb0s = []
    p_cb0s = []
    mean_isis_buffer = []
    cv_isis_buffer = []
    for bT in [10]:
        cb0 = (1 * 0.2 * bT) / (1 * 0.2 + 5)
        p_cb0 = cb0 / (0.2 + cb0)
        print(p_cb0)
        cb0s.append(cb0)
        p_cb0s.append(p_cb0)
        folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin_buffer/out/"
        file_spikes = f"spike_times_langevin_buffer_bt{bT:.2f}_ip1.00_tau5.62e+00_j1.26e-02_K10_5.dat"
        file_dynamics = f"ca_langevin_bt{bT:.2f}_ip1.00_tau5.62e+00_j1.26e-02_K10_5.dat"
        isis = np.loadtxt(folder + file_spikes)
        mean_isi = np.mean(isis)
        cv_isi = np.std(isis)/mean_isi
        data = np.loadtxt(folder + file_dynamics)
        ts, cas, _, _, _, cbs = np.transpose(data)

        ax1.set_ylabel(r"$c_{\rm i}$")
        ax1.set_xlabel("$t$ /s")
        ax1.plot(ts, cas, lw=1, color=st.colors[0], label=rf"$\langle T \rangle = {mean_isi:.0f}, CV_T = {cv_isi:.1f}$")
        ax1.set_xlim([0, 1000])
        ax1.legend(fancybox=False, fontsize=9)

        p_cbs = []
        for i, (ca, cb) in enumerate(zip(cas, cbs)):
            p_cb = cb/(cb + ca)
            p_cbs.append(p_cb)
        ax2.set_xlabel("$t$ /s")
        ax2.set_ylabel(r"$c_{b}/c_T$")
        ax2.plot(ts, cbs, lw=1, color=st.colors[4])
        ax2.set_xlim([0, 1000])

    for bT in [50]:
        cb0 = (1 * 0.2 * bT) / (1 * 0.2 + 5)
        p_cb0 = cb0 / (0.2 + cb0)
        print(p_cb0)
        cb0s.append(cb0)
        p_cb0s.append(p_cb0)
        folder = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin_buffer/out/"
        file_spikes = f"spike_times_langevin_buffer_bt{bT:.2f}_ip1.00_tau5.62e+00_j1.26e-02_K10_5.dat"
        file_dynamics = f"ca_langevin_bt{bT:.2f}_ip1.00_tau5.62e+00_j1.26e-02_K10_5.dat"
        isis = np.loadtxt(folder + file_spikes)
        mean_isi = np.mean(isis)
        cv_isi = np.std(isis)/mean_isi
        data = np.loadtxt(folder + file_dynamics)
        ts, cas, _, _, _, cbs = np.transpose(data)

        ax3.set_ylabel(r"$c_{\rm i}$")
        ax3.set_xlabel("$t$ /s")
        ax3.plot(ts, cas, lw=1, color=st.colors[0], label=rf"$\langle T \rangle = {mean_isi:.0f}, CV_T = {cv_isi:.1f}$")
        ax3.set_xlim([0, 1000])
        ax3.legend(fancybox=False, fontsize=9)

        p_cbs = []
        for i, (ca, cb) in enumerate(zip(cas, cbs)):
            p_cb = cb/(cb + ca)
            p_cbs.append(p_cb)
        ax4.set_xlabel("$t$ /s")
        ax4.set_ylabel(r"$c_{b}/c_T$")
        ax4.plot(ts, cbs, lw=1, color=st.colors[4])
        ax4.set_xlim([0, 1000])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig9b.pdf", transparent=True)

    plt.show()