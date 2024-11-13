import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as stats
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    # Set up plot style
    st.set_default_plot_style()
    w = 6.5*1.25
    h = 3.*1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    axm1 = fig.add_subplot(gs[1, 2])
    axm2 = fig.add_subplot(gs[0, 2])
    gs11 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0, 0], wspace=0.1)
    gs12 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0, 1], wspace=0.1)
    gs21 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, 0], wspace=0.1)
    gs22 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, 1], wspace=0.1)
    ax1 = fig.add_subplot(gs11[0:2])
    ax2 = fig.add_subplot(gs11[2])
    ax3 = fig.add_subplot(gs12[0:2])
    ax4 = fig.add_subplot(gs12[2])
    ax5 = fig.add_subplot(gs21[0:2])
    ax6 = fig.add_subplot(gs21[2])
    ax7 = fig.add_subplot(gs22[0:2])
    ax8 = fig.add_subplot(gs22[2])
    st.remove_top_right_axis([axm1, axm2, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8])

    axm1.text(0.1, 0.95, r"D", fontsize=11, transform=axm1.transAxes, va='top')
    axm2.text(0.1, 0.95, r"C", fontsize=11, transform=axm2.transAxes, va='top')
    ax1.text(0.1, 0.95, r"A$_1$", fontsize=11, transform=ax1.transAxes, va='top', bbox=dict(facecolor='w', edgecolor='none'))
    ax3.text(0.1, 0.95, r"A$_2$", fontsize=11, transform=ax3.transAxes, va='top', bbox=dict(facecolor='w', edgecolor='none'))
    ax5.text(0.1, 0.95, r"B$_1$", fontsize=11, transform=ax5.transAxes, va='top', bbox=dict(facecolor='w', edgecolor='none'))
    ax7.text(0.1, 0.95, r"B$_2$", fontsize=11, transform=ax7.transAxes, va='top', bbox=dict(facecolor='w', edgecolor='none'))

    for ax in [ax1, ax5]:
        ax.set_ylabel(rf"$Y(t; \Delta t = 0.1)$")
    for ax in [ax3, ax7]:
        ax.set_ylabel(rf"$Y(t; \Delta t = 1.0)$")
    for ax in [ax2, ax4, ax6, ax8]:
        ax.set_xlabel("$P(Y)$")
        ax.set_yticks([])
        ax.set_xticks([])
    for ax in [ax1, ax3, ax5, ax7]:
        ax.set_xlabel("$t$ / s")
        ax.set_xlim([0, 10])
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_ylim([0, 5])
    for ax in [ax5, ax6, ax7, ax8]:
        ax.set_ylim([0, 20])

    for ax in [axm1, axm2]:
        ax.set_xlabel("$\Delta t$ / s")
        ax.set_xscale("log")
    axm1.set_yscale("log")
    axm2.set_ylabel("$\gamma_Y$")
    axm1.set_ylabel("$t$ / s")

    color = st.Colors

    axm1.set_ylim([0.01, 1000])
    axm1.set_yticks([0.01, 0.1, 1, 10, 100, 1000])
    axm1.set_xlim([0.01, 1])
    axm1.plot([0.01, 1], [0.01, 1], ls=":", color="k")
    axm1.text(0.1, np.sqrt(10), r"$\tau$", fontsize=13, va="center", ha="center")
    axm1.text(0.1, 10 * np.sqrt(10), r"$\langle T \rangle$", fontsize=13, va='center', ha="center")
    axm1.text(0.1, 0.04, r"$\tau_Y$", fontsize=13, va="center", ha="center")
    axm1.fill_between([0.01, 1], [1, 1], [10, 10], hatch="///", edgecolor=color.green[1], facecolor="none")
    axm1.fill_between([0.01, 1], [10, 10], [100, 100], hatch="///", edgecolor=color.green[1], facecolor="none")
    #axm1.arrow(0.02, 0.02, 0, 1 - 0.02, fc="k", lw=0.25, zorder=1)

    # Set up Parameters
    ca1 = 0.20
    ca2 = 0.50
    for i, ca in enumerate([ca1, ca2]):
        if i == 0:
            plt_color = color.palette[0]
        else:
            plt_color = color.palette[2]
            
        ca_rest = 0.2
        K = 10
        N = 5
        M = 3
        jca = 1
        tau = 1
        # Load Data
        home = os.path.expanduser("~")
        folder = "/Data/calcium_spikes_markov/ca_fix/"
        file = f"ca_markov_cafix{ca:.2f}_ip1.00_tau{tau:.2e}_j{jca:.2e}_K10_5.dat"
        data = np.loadtxt(home + folder + file)
        ts, cas, jpuffs, adaps = np.transpose(data)
        dt = ts[1] - ts[0]

        file_puff = f"puff_markov_cafix{ca:.2f}_ip1.00_tau{tau:.2e}_j{jca:.2e}_K1_5.dat"
        data_puff = np.loadtxt(home + folder + file_puff)
        data_puff_spaced = fc.equally_spaced_puff_data(data_puff, dt=0.01)
        ns_puff_spaced = [x[1] for x in data_puff_spaced]


        corr_taus = []
        corr_times = []
        corr_func = fc.autocorrelation_function(ns_puff_spaced, dt=0.01)
        corr_time = (sum(corr_func[1]) * 0.01) / np.var(ns_puff_spaced)
        corr_taus.append(0.01)
        corr_times.append(corr_time)
        for f in [2, 3, 4, 5, 6, 8, 9, 10, 16, 32, 64, 100]:
            print(f)
            ns_puff_coarse_grained = fc.moving_coarse_grain_list(ns_puff_spaced, f)
            corr_func = fc.autocorrelation_function(ns_puff_coarse_grained, dt=0.01*f)
            corr_time = (sum(corr_func[1]) * 0.01) / np.var(ns_puff_coarse_grained)
            corr_taus.append(0.01*f)
            corr_times.append(corr_time)
        print(corr_times)
        if ca==ca1:
            axm1.plot(corr_taus, corr_times, c=plt_color, zorder=3)
        else:
            axm1.plot(corr_taus, corr_times, c=plt_color, zorder=2)

        r_opn_s = 0.1
        r_ref = 20
        r_cls = 50

        delta_t = ts[1] - ts[0]
        print(delta_t)
        dt1 = 0.1
        dt2 = 1.0
        n1 = int(dt1 / (2 * delta_t))
        n2 = int(dt2 / (2 * delta_t))

        jpuffs_f1 = fc.moving_coarse_grain_list(jpuffs, n1)
        jpuffs_f2 = fc.moving_coarse_grain_list(jpuffs, n2)
        moment3s = []
        std_moment3s = []
        fs = np.arange(1, 100, 2)
        for f in fs:
            print(f)
            jpuffs_f = fc.moving_coarse_grain_list(jpuffs, f)
            ms = []
            for k in range(10):
                l = int(len(jpuffs_f)/10)
                ms.append(stats.skew(jpuffs_f[k*l:(k+1)*l]))
            moment3_f = np.mean(ms)
            std_moment3_f = np.std(ms)
            moment3s.append(moment3_f)
            std_moment3s.append(std_moment3_f)

        mean_f1 = np.mean(jpuffs_f1)
        mean_f2 = np.mean(jpuffs_f2)
        std_f1 = np.std(jpuffs_f1)
        std_f2 = np.std(jpuffs_f2)
        covar_f1 = fc.k_corr(jpuffs_f1, jpuffs_f1, k=1)
        covar_f2 = fc.k_corr(jpuffs_f2, jpuffs_f2, k=1)
        print(covar_f1, covar_f2)
        gaus1 = np.linspace(mean_f1 - 3 * std_f1, mean_f1 + 3 * std_f1, 100)
        gaus2 = np.linspace(mean_f2 - 3 * std_f2, mean_f2 + 3 * std_f2, 100)

        axm2.set_ylim([0, 3])
        axm2.set_xlim([0.01, 1])
        if ca == ca1:
            label=r"$c_R$"
        else:
            label=r"$c_T$"
        axm2.plot([delta_t*f for f in fs], moment3s, c=plt_color, label=label)
        #axm2.fill_between([delta_t*f for f in fs], [x + s/np.sqrt(10) for (x, s) in zip(moment3s, std_moment3s)], [x - s/np.sqrt(10) for (x, s) in zip(moment3s, std_moment3s)], color=color.palette[i], alpha=0.5)

        if ca == ca1:
            mean_jpuff = K*fc.mean_jp_single_theory(ca1, N, M, 1)
            d_jpuff = K*fc.noise_intensity_jp_single_theory(ca1, N, M, 1)
            std_jpuff_f1 = np.sqrt(2*d_jpuff/dt1)
            std_jpuff_f2 = np.sqrt(2*d_jpuff/dt2)
            if n1==0:
                ax1.plot(ts, jpuffs_f1, c=plt_color)
            else:
                ax1.plot(ts[n1:-n1], jpuffs_f1, c=plt_color)
            #ax1.axhline(mean_jpuff, ls="--", color=color.palette[5])

            ax3.plot(ts[n2:-n2], jpuffs_f2, c=plt_color)
            #ax3.axhline(mean_jpuff, ls="--", color=color.palette[5])

            ax2.plot(fc.gaussian_dist_std(gaus1, mean_jpuff, std_jpuff_f1), gaus1, c=color.palette[5], lw=1.)
            ax2.hist(jpuffs_f1, bins=25, alpha=0.5, color=plt_color, density=True, orientation="horizontal")

            ax4.plot(fc.gaussian_dist_std(gaus2, mean_jpuff, std_jpuff_f2), gaus2, c=color.palette[5], lw=1.)
            ax4.hist(jpuffs_f2, bins=50, alpha=0.5, color=plt_color, density=True, orientation="horizontal")

        if ca == ca2:
            mean_jpuff = K*fc.mean_jp_single_theory(ca2, N, M, 1)
            print("Calculate intensity...")
            d_jpuff = K*fc.noise_intensity_jp_single_theory(ca2, N, M, 1)
            print("...done")
            std_jpuff_f1 = np.sqrt(2*d_jpuff/dt1)
            std_jpuff_f2 = np.sqrt(2*d_jpuff/dt2)
            print(std_jpuff_f1, np.std(jpuffs_f1))
            print(std_jpuff_f2, np.std(jpuffs_f2))
            print()
            if n1==0:
                ax5.plot(ts, jpuffs_f1, c=plt_color)
            else:
                ax5.plot(ts[n1:-n1], jpuffs_f1, c=plt_color)
            #ax5.axhline(mean_jpuff, ls="--", color=color.palette[5])

            ax7.plot(ts[n2:-n2], jpuffs_f2, c=plt_color)
            #ax7.axhline(mean_jpuff, ls="--", color=color.palette[5])

            ax6.plot(fc.gaussian_dist_std(gaus1, mean_jpuff, std_jpuff_f1), gaus1, c=color.palette[5], lw=1.)
            ax6.hist(jpuffs_f1, bins=25, alpha=0.5, color=plt_color, density=True, orientation="horizontal")

            ax8.plot(fc.gaussian_dist_std(gaus2, mean_jpuff, std_jpuff_f2), gaus2, c=color.palette[5], lw=1.)
            ax8.hist(jpuffs_f2, bins=25, alpha=0.5, color=plt_color, density=True, orientation="horizontal")
    axm2.legend(fancybox=False)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/SUB2/figures/fig4.png", transparent=True, dpi=300)
    plt.show()
