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
    taus = [5, 1]
    ps = [0.015, 0.060]
    tau_ers = np.logspace(1, 3, 21)
    tau_er = tau_ers[13]
    print(tau_er)
    eps_er = 0.1

    cmap_YlGnBu = plt.get_cmap("YlGnBu", 10)
    for tau, p, ax, axr in zip(taus, ps, [ax1a, ax2a], [ax1b, ax2b]):
        data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)

        mean = np.mean(data_isi)
        std = np.std(data_isi)
        cv = std/mean
        cv2 = cv**2
        ts = np.linspace(mean-3*std, mean+3*std, 501)
        p_inv_gaus = fc.inverse_gaussian_dist(ts, mean, cv2)
        p_gamma = fc.gamma_dist(ts, mean, cv2)
        l1, = ax.plot(ts, p_inv_gaus, lw=1, c="k", ls="--", label="Inv.\ Gaus.")
        l2, = ax.plot(ts, p_gamma, lw=1, c="k", ls=":", label="Gamma")
        ax.hist(data_isi, bins=50, color=st.colors[1], alpha=0.75, density=True)
        legend1 = ax.legend(loc=1, handles=[l1, l2], fancybox=False, framealpha=1., edgecolor="k", fontsize=8)
        ax.add_artist(legend1)

        c1_patch = mpatches.Patch(color=st.colors[1], alpha=0.75, label=rf"$\langle T \rangle = {mean:.0f}$" + "\n" + rf"$CV_T={cv:.1f}$")
        legend2 = ax.legend(loc=5, handles=[c1_patch], fancybox=False, framealpha=1., edgecolor="k", fontsize=8)
        ax.add_artist(legend2)

        p1 = fc.k_corr(data_isi, data_isi, 1)/(std**2)
        h1 = axr.hist2d(data_isi[:-1], data_isi[1:], density=True, bins=25, cmap=cmap_YlGnBu)
        axr.set_xlabel(r"$T_i$")
        axr.set_ylabel(r"$T_{i+1}$")
        c2_patch = mpatches.Patch(color=st.colors[1], alpha=0.75, label=rf"$\rho_1 = {p1:.2f}$")
        legend3 = axr.legend(loc=1, handles=[c2_patch], fancybox=False, framealpha=1., edgecolor="k", fontsize=8)

        divider = make_axes_locatable(axr)
        caxr = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(h1[3], cax=caxr, orientation='vertical')
        cbar.set_label(r"$p(T_i, T_{i+1})$", loc="center")
    plt.show()