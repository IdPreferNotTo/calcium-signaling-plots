import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    n = 16
    print(np.logspace(1, 3, 50)[-n])
    taua = 244
    ampa = 0.11
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap"
    file_isi = f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    isis = np.loadtxt(folder + file_isi)
    isis = [I for I in isis[100:]]
    mean = np.mean(isis)
    std = np.std(isis)
    var = np.var(isis)
    cv = std/mean
    print(mean)

    ks = np.arange(1, 5)
    rho_ks = []
    for k in ks:
        var_k_isis = fc.k_corr(isis, isis, k)
        rho_k = var_k_isis/var
        rho_ks.append(rho_k)

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
    ax0 = fig.add_subplot(gs1[1:3, 0:2])
    ax1 = fig.add_subplot(gs1[1:3, 2])
    ax2 = fig.add_subplot(gs1[0, 0:2])
    ax3 = fig.add_subplot(gs2[0])
    ax4 = fig.add_subplot(gs2[1:3])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3, ax4])

    ax0.set_xlabel("$T_i$")
    ax0.set_ylabel("$T_{i+1}$")
    ax0.set_xlim([0.5 * mean, 1.5 * mean])
    ax0.set_ylim([0.5 * mean, 1.5 * mean])
    ax0.scatter(isis[:-1], isis[1:], fc="w", ec=st.colors[1], s=20, zorder=3)
    xmin, xmax = ax0.get_xlim()
    ymin, ymax = ax0.get_ylim()
    xs = np.linspace(0.6*mean, 1.4*mean)
    ys = mean + rho_ks[0]*(xs - mean)
    ax0.plot(xs, ys, c="k", zorder=5)

    ax1.set_xlabel("$P_0(T_{i+1})$")
    ax1.set_ylim([ymin, ymax])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.hist(isis, bins=20, color=st.colors[1], density=True, alpha=0.6, orientation="horizontal")
    ts_inv_gaussian, ps_inv_gaussian = fc.inverse_gaussian(mean, cv)
    ax1.plot(ps_inv_gaussian, ts_inv_gaussian, color="k", ls="--")

    ax2.set_ylabel("$P_0(T_{i})$")
    ax2.set_xlim([xmin, xmax])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.hist(isis, bins=20, color=st.colors[1], density=True, alpha=0.6)
    ax2.plot(ts_inv_gaussian, ps_inv_gaussian, color="k", ls="--")

    ax3.set_ylim([-0.3, 0.05])
    ax3.set_xlabel("$k$")
    ax3.set_ylabel(r"$\rho_k$")
    ax3.scatter(ks, rho_ks, fc="w", ec=st.colors[1], s=20, zorder=3)
    ax3.axhline(0, ls=":", c="C7")

    ax4.set_xlabel("$\omega$")
    ax4.set_ylabel("$S(\omega)$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_xlim([0.001, 2])
    ax4.plot(ws, spectrum, color=st.colors[1], label="unshuffled")
    ax4.plot(ws, spectrum_shuffle, color=st.colors[0], label="shuffled")
    ax4.legend(frameon=False, loc=4)

    ax4.axhline((1. / mean), ls="--", c="C7")
    ax4.axhline((1. / mean) * cv ** 2, ls="--", c="C7")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig4.png", transparent=True)
    plt.show()