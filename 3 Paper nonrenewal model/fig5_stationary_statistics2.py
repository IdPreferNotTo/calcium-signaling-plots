import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

import functions as fc
import styles as st
import default_parameters as df


def conditional_mean(T1s, Tmin, Tmax, num):
    dT = (Tmax - Tmin)/num
    T2s = [[] for _ in range(num)]
    for T1, T2 in zip(T1s[:-1], T1s[1:]):
        idx = (T1 - (Tmin - dT/2))/dT
        if idx < 0 or idx >= num:
            continue
        T2s[int(idx)].append(T2)
    return np.linspace(Tmin, Tmax, num), [np.mean(T2) for T2 in T2s]


if __name__ == "__main__":
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 1.55 * 1.25
    # fig = plt.figure(tight_layout=True, figsize=(w, h))
    # gs = gridspec.GridSpec(1, 2)
    # ax1 = fig.add_subplot(gs[0])
    # ax2 = fig.add_subplot(gs[1])
    fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(w, h))
    ax1 = axs[0]
    ax2 = axs[1]
    axis = [ax2]
    st.remove_top_right_axis(axis)
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')

    ax1.set_title(r"$p(T_{i+1}, T_{i})$")
    ax1.set_xlabel(r"$T_i$ / s")
    ax1.set_ylabel(r"$T_{i+1}$ / s")
    # ax1.set_yticks([1, 2])
    # ax1.set_xticks([1, 2])

    axins1 = ax1.inset_axes((0.55, 0.5, .4, .4))
    axins1.set_title(r"$\rho_k$", fontsize=8)
    axins1.set_xlabel(r"$k$", fontsize=8)
    axins1.set_xticks([1, 2, 3, 4, 5])
    axins1.set_ylim([-0.6, 0.15])
    axins1.axhline(ls=":", lw=1, c="C7")
    axins1.tick_params(which="both", direction="in", labelsize=7)

    ax2.set_xlabel("$f$ / Hz")
    ax2.set_ylabel("$S(f)$ / Hz")
    ax2.set_ylim([0, 0.032])
    ax2.set_xlim([0, 0.032])
    axins2 = ax2.inset_axes((0.55, 0.5, .40, .40))

    axins2.tick_params(which="both", direction="in")
    axins2.set_title("$S(f)$", fontsize=8)
    axins2.set_xscale("log")
    axins2.set_yscale("log")
    axins2.set_xlim([0.0001, 0.005])
    axins2.set_ylim([0.0001, 0.005])
    axins2.set_yticklabels([])
    axins2.set_xticklabels([])

    rect = patches.Rectangle((0., 0.), 0.005, 0.005, linewidth=0.75, edgecolor='k', facecolor='none', zorder=1)
    ax2.add_patch(rect)
    ax2.plot([0.0, 0.55], [0.005/0.032, 0.90], ls=":", c="k", lw=0.75, transform=ax2.transAxes)
    ax2.plot([0.005/0.032, 0.955], [0.0, 0.50], ls=":", c="k", lw=0.75, transform=ax2.transAxes)
    # Parameters
    tau = 5
    p = 0.015
    tau_er = 300
    eps_er = 0.03

    cmap_YlGnBu = plt.get_cmap("YlGnBu", 10)

    home = os.path.expanduser("~")
    file = home + f"/Data/calcium/markov/adaptive/sequence/long_spike_times_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
    data_isi = np.loadtxt(file)
    print(len(data_isi))
    mean = np.mean(data_isi)
    std = np.std(data_isi)

    h1 = ax1.hist2d(data_isi[:-1], data_isi[1:], range=[[0.5*mean, 2.1*mean], [0.5*mean, 2.1*mean]], density=True, bins=25, cmap=cmap_YlGnBu)
    cbar = fig.colorbar(h1[3], ax=ax1, orientation='vertical')
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.yaxis.set_offset_position('left')

    n_min = 0.7*mean
    n_max = 1.3*mean
    dn = (n_max - n_min) / 20
    norm_isi = data_isi[1:]
    con_isi = [[] for _ in range(20)]
    for n_isi2, n_isi1 in zip(norm_isi[1:], norm_isi[:-1]):
        idx = int((n_isi1 - n_min)/dn)
        if idx < 0 or idx > 19:
            continue
        else:
            con_isi[idx].append(n_isi2)
    con_mean_isi = [np.mean(isis) for isis in con_isi]
    ax1.plot(np.linspace(n_min, n_max, 20), con_mean_isi, c="r")
    ax1.axhline(mean, c="k", ls=":")
    ax1.axvline(mean, c="k", ls=":")

    Tmax = 5000
    fmin = 1/Tmax
    fmax = 0.035
    ks = np.arange(1, 6)
    pks = []
    for k in ks:
        pks.append(fc.k_corr(data_isi, data_isi, k)/(std**2))
    axins1.scatter(ks, pks, fc="w", ec=st.colors[0], s=20, zorder=3)
    axins1.plot(ks, pks[0]*np.exp(-(ks-1)/1.346), c=st.colors[0])
    print(pks[0])

    n = int(fmax/fmin)
    fs = fmin*np.arange(n)
    Ps = fc.power_spectrum_isis(fs, data_isi, Tmax=Tmax)
    ax2.plot(fs, Ps, c=st.colors[1], zorder=3)
    axins2.plot(fs, Ps, c=st.colors[1])

    np.random.shuffle(data_isi)
    Tmax = 5000
    fmin = 1/Tmax
    fmax = 0.035
    ks = np.arange(1, 6)
    pks = []
    for k in ks:
        pks.append(fc.k_corr(data_isi, data_isi, k)/(std**2))
    axins1.scatter(ks, pks, fc="w", ec=st.colors[5], s=20, zorder=2)
    n = int(fmax/fmin)
    fs = fmin*np.arange(n)
    Ps2 = fc.power_spectrum_isis(fs, data_isi, Tmax=Tmax)
    ax2.plot(fs, Ps2, c=st.colors[5], zorder=2)
    axins2.plot(fs, Ps2, c=st.colors[5])


    home = os.path.expanduser("~")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig5.pdf", dpi=300, transparent=True)
    plt.show()