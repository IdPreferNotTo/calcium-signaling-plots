import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    taua = 829
    ampa = 0.0309
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap"
    file_isi = f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    isis = np.loadtxt(folder + file_isi)
    isis = [I for I in isis[100:]]
    mean_isis = np.mean(isis)
    std_isis = np.std(isis)
    var_isis = np.var(isis)

    ks = np.arange(1, 5)
    rho_ks = []
    for k in ks:
        var_k_isis = fc.k_corr(isis, isis, k)
        rho_k = var_k_isis/var_isis
        rho_ks.append(rho_k)

    mean = np.mean(isis)
    cv = np.std(isis)/mean

    ws = np.logspace(-3, 0, 100)
    spectrum = fc.power_spectrum_isis(ws, isis, 2000)
    isis_shuffled = list(isis)
    np.random.shuffle(isis_shuffled)
    spectrum_shuffle = fc.power_spectrum_isis(ws, isis_shuffled, 2000)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3*2.25, 3))
    gs = gridspec.GridSpec(1, 2)
    gs1 = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs[0], hspace=0.0, wspace=0.0)
    gs2 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1], hspace=1.)
    ax0 = fig.add_subplot(gs1[1:3, 2])
    ax1 = fig.add_subplot(gs1[1:3, 0:2])
    ax2 = fig.add_subplot(gs1[0, 0:2])
    ax3 = fig.add_subplot(gs2[0])
    ax4 = fig.add_subplot(gs2[1:3])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3, ax4])

    ax0.set_xlabel("$P_0(T_{i+1})$")
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax0.hist(isis, bins=20, color=st.colors[1], density=True, alpha=0.6, orientation="horizontal")

    ax1.set_xlabel("$T_i$")
    ax1.set_ylabel("$T_{i+1}$")
    ax1.scatter(isis[:-1], isis[1:], fc="w", ec=st.colors[1], s=20, zorder=3)

    ax2.set_ylabel("$P_0(T_{i})$")
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.hist(isis, bins=20, color=st.colors[1], density=True, alpha=0.6)

    ax3.set_ylim([-0.1, 0.1])
    ax3.set_xlabel("$k$")
    ax3.set_ylabel(r"$\rho_k$")
    ax3.scatter(ks, rho_ks, fc="w", ec=st.colors[1], s=20, zorder=3)

    ax4.set_xlabel("$\omega$")
    ax4.set_ylabel("$S(\omega)$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.plot(ws, spectrum, color=st.colors[1])
    ax4.plot(ws, spectrum_shuffle, color=st.colors[0])

    ax4.axhline((1./mean), ls="--", c="C7")
    ax4.axhline((1./mean)*cv ** 2 , ls="--", c="C7")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig3.pdf", transparent=True)
    plt.show()