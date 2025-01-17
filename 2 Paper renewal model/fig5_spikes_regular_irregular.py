import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import scipy.special as special
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    w = 6.5*1.25
    h = 3.*1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    gs11 = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs[0, 0], hspace=0.1)
    gs12 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0, 1], wspace=0.1)
    gs13 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs[0, 2], wspace=0.05)

    gs21 = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs[1, 0], hspace=0.1)
    gs22 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[1, 1], wspace=0.1)
    gs23 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs[1, 2], hspace=0.05)

    ax11= fig.add_subplot(gs11[0:2])
    ax12 = fig.add_subplot(gs11[2])
    ax2 = fig.add_subplot(gs12[0])
    axin2 = ax2.inset_axes([0.65, 0.65, 0.35, 0.35])
    ax3 = fig.add_subplot(gs13[:])
    ax3a = fig.add_subplot(gs13[0])
    ax3b = fig.add_subplot(gs13[1])#ax3.inset_axes([0.45, 0.45, 0.5, 0.5])

    ax41 = fig.add_subplot(gs21[0:2])
    ax42 = fig.add_subplot(gs21[2])
    ax5 = fig.add_subplot(gs22[0])
    axin5 = ax5.inset_axes([0.65, 0.65, 0.35, 0.35])
    ax6  = fig.add_subplot(gs23[:])
    ax6a = fig.add_subplot(gs23[0])
    ax6b = fig.add_subplot(gs23[1]) #ax6.inset_axes([0.45, 0.45, 0.5, 0.5])

    st.remove_everything([ax3, ax6])
    st.remove_top_right_axis([ax11, ax12, ax2, axin2, ax3a, ax3b, ax41, ax42, ax5, axin5, ax6a, ax6b])
    axis1 = [ax11, ax12, ax2, axin2, ax3a, ax3b]
    axis2 = [ax41, ax42, ax5, axin5, ax6a, ax6b]
    axiss = [axis1, axis2]

    ax11.text(0.05, 0.95, r"A$_1$", fontsize=11, transform=ax11.transAxes, va='top')
    ax2.text(0.05, 0.95, r"A$_2$", fontsize=11, transform=ax2.transAxes, va='top')
    ax3a.text(0.10, 0.95, r"A$_3$", fontsize=11, transform=ax3a.transAxes, va='top')
    ax41.text(0.05, 0.95, r"B$_1$", fontsize=11, transform=ax41.transAxes, va='top')
    ax5.text(0.05, 0.95, r"B$_2$", fontsize=11, transform=ax5.transAxes, va='top')
    ax6a.text(0.05, 0.95, r"B$_3$", fontsize=11, transform=ax6a.transAxes, va='top')
    home = os.path.expanduser("~")
    # Plot regular spiketrain
    taus = [5, 1] #[5.62, 1.78]
    js = [0.015, 0.06] #[0.0126, 0.0355]
    for k, (tau, j, axis) in enumerate(zip(taus, js, axiss)):

        ax11, ax12, ax2, axin2, ax3a, ax3b = axis
        ca_r = 0.20
        ca_t = 0.50
        N = 10
        n = 5
        m = 3
        ip3 = 1.0
        r_opn_single = 0.1
        r_ref = 20

        # Get Data
        folder = home + "/Data/calcium_spikes_markov/Data_no_adap/"
        file_calcium = f"ca_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        data_calcium = np.loadtxt(folder + file_calcium)
        data_isis = np.loadtxt(folder + file_spikes)
        print(len(data_isis))
        ts, cas, jpuffs, adaps = np.transpose(data_calcium)

        folder_puffs = home + "/Data/calcium_spikes_markov/Data_no_adap/"
        file_puffs = f"puff_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        data_puff = np.loadtxt(folder_puffs + file_puffs)
        ts_puff, _, _, ns_puff = np.transpose(data_puff)

        file_p0_theory = home + f"/Data/calcium_spikes_theory/probability_density/p0_ip1.00_tau{tau:.2e}_j{j:.2e}_K10.dat"
        data = np.loadtxt(file_p0_theory)
        ca_theory, p0_theory = np.transpose(data)

        file_r0_theory = home + f"/Data/calcium_spikes_theory/firing_rate/r0_ip1.00_tau{tau:.2e}_j{j:.2e}_K10.dat"
        r0 = np.loadtxt(file_r0_theory)
        mean_isi_langevin = 1/r0

        file_cv_theory = home + f"/Data/calcium_spikes_theory/coefficient_of_variation/cv_ip1.00_tau{tau:.2e}_j{j:.2e}_K10.dat"
        cv_isi_langevin = np.loadtxt(file_cv_theory)


        if k == 0:
            plt_color = color.palette[0]
        else:
            plt_color = color.palette[2]
            
        #Get drift
        xs = np.linspace(0.2, 0.5, 100)
        dx = xs[1] - xs[0]
        ca_drifts = []
        for x in xs:
            ca_drift = -(x- ca_r)/tau + j*N*fc.mean_jp_single_theory(x, n, m, ip3)
            if ca_drift < 0 and ca_drift_tmp > 0:
                ca_fix = x
            ca_drift_tmp = ca_drift
            ca_drifts.append(20*ca_drift)

        # Plots ax11 ax12
        count = 0
        spike_times = []
        max_i = 0
        max_spikes = 5
        while count < max_spikes:
            if cas[max_i] == 0.5:
                count += 1
                spike_times.append(ts[max_i])
            max_i += 1

        t_offset = ts[500]
        ts_shift = ts[500:] - ts[500]
        ts_puff_shift = ts_puff[500:] - ts_puff[500]
        cas_shift = cas[500:]

        ax11.plot(ts_shift, cas_shift, lw=1, color=plt_color)
        ax12.plot(ts_puff_shift, [j*n for n in ns_puff[500:]], lw=1, color=plt_color)
        #ax12.plot(ts[:max_i], jpuffs[:max_i], color=plt_color)

        #Plot ax2
        cas = [ca for ca in cas if ca !=1]
        ax2.hist(cas, bins=25, density=True, alpha=0.5, color=plt_color)
        ax2.plot(ca_theory[::100], p0_theory[::100], color=color.palette[5])
        axin2.plot(xs, ca_drifts, color=plt_color)
        axin2.set_title(r"$f(c_{\rm i})$", fontsize=11)
        axin2.set_yticks([0])
        axin2.axhline(0, lw=1, ls=":", color="C7")
        if k == 0:
            axin2.set_ylim(bottom=-0.05)

        #Plot ax3
        mean_isi = np.mean(data_isis)
        cv_isi = np.std(data_isis) / mean_isi
        print(mean_isi, cv_isi)
        #print(mean_isi, mean_isi_theory, cv_isi, cv_theory)
        if k == 0:
            ax3a.hist(data_isis, bins=20, color=plt_color, density=True, alpha=0.5)
            ax3b.hist(data_isis, bins=20, color=plt_color, density=True, alpha=0.5)
            ax3a.text(0.7, 0.7, r"inv.\ Gaus.", fontsize=9, color=color.palette[5], transform=ax3a.transAxes,
                      va='center', ha='center')
            ax3b.text(0.7, 0.7, r"Gamma", fontsize=9, color=color.palette[5], transform=ax3b.transAxes, va='center',
                      ha='center')
        else:
            ax3a.hist(data_isis, bins=50, color=plt_color, density=True, alpha=0.5)
            ax3b.hist(data_isis, bins=50, color=plt_color, density=True, alpha=0.5)
            ax3a.text(0.7, 0.5, r"inv.\ Gaus.", fontsize=9, color=color.palette[5], transform=ax3a.transAxes,
                      va='center', ha='center')
            ax3b.text(0.7, 0.5, r"Gamma", fontsize=9, color=color.palette[5], transform=ax3b.transAxes, va='center',
                      ha='center')

        ts_dist = np.linspace(1, 150, 1001)
        inv_gaus_dist_markov = []
        inv_gaus_dist_langevin =[]
        for t in ts_dist:
            p_markov = np.sqrt(mean_isi / (2 * np.pi * np.power(cv_isi, 2) * (t ** 3))) * np.exp(
                -(t - mean_isi) ** 2 / (2 * mean_isi * np.power(cv_isi, 2) * t))
            p_langevin = np.sqrt(mean_isi_langevin / (2 * np.pi * np.power(cv_isi_langevin, 2) * (t ** 3))) * np.exp(
                -(t - mean_isi_langevin) ** 2 / (2 * mean_isi_langevin * np.power(cv_isi_langevin, 2) * t))
            inv_gaus_dist_markov.append(p_markov)
            inv_gaus_dist_langevin.append(p_langevin)

        ax3a.plot(ts_dist, inv_gaus_dist_markov, color=color.palette[5], ls="--", label="Inv.\ Gaussian")
        ax3a.plot(ts_dist, inv_gaus_dist_langevin, color=color.palette[5], ls=":")


        ts_dist = np.linspace(1, 150, 1001)
        gamma_dist_markov = []
        gamma_dist_langevin =[]
        alpha_markov = 1/np.power(cv_isi,2)
        beta_markov  = 1/(mean_isi*np.power(cv_isi,2))
        alpha_langevin = 1/np.power(cv_isi_langevin, 2)
        beta_langevin = 1/(mean_isi_langevin*np.power(cv_isi_langevin,2))

        for t in ts_dist:
            tscale = t / (mean_isi * np.power(cv_isi, 2))
            tsclae_langevin = t /(mean_isi_langevin * np.power(cv_isi_langevin, 2))
            p_markov = np.power(tscale, 1/cv_isi**2) * np.exp(-tscale)/(t * special.gamma(1/cv_isi**2))
            p_langevin = np.power(tsclae_langevin, 1/cv_isi_langevin**2) * np.exp(-tsclae_langevin)/(t * special.gamma(1/cv_isi_langevin**2))
            gamma_dist_markov.append(p_markov)
            gamma_dist_langevin.append(p_langevin)
        #gamma_dist_markov = stats.gamma.pdf(ts_dist, a=alpha_markov, scale=1 / beta_markov)
        #gamma_dist_langevin = stats.gamma.pdf(ts_dist, a = alpha_langevin, scale = 1/ beta_langevin)

        ax3b.plot(ts_dist, gamma_dist_markov, color=color.palette[5], ls="--", label="Gamma")
        ax3b.plot(ts_dist, gamma_dist_langevin, color=color.palette[5], ls=":")

        for i in range(max_spikes-1):
            x_left = spike_times[i] - t_offset
            x_right = spike_times[i+1] -t_offset
            dx = spike_times[i+1] - spike_times[i]
            if i == 1:
                ax11.text(x_left + dx/2, 0.65, "$T_{i-1}$", ha="center", va="center", clip_on=False)
            elif i == 2:
                ax11.text(x_left + dx/2, 0.65, "$T_{i}$", ha="center", va="center", clip_on=False)
            elif i == 3:
                ax11.text(x_left + dx / 2, 0.65, "$T_{i+1}$", ha="center", va="center", clip_on=False)
            else:
                ax11.text(x_left + dx / 2, 0.65, "$\dots$", ha="center", va="center", clip_on=False)
            ax11.arrow(x_left + 0.05*dx, 0.55, 0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)
            ax11.arrow(x_right -0.05*dx, 0.55, -0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)

        #ax0.set_xlabel("$t$ / s")
        ax11.set_ylabel(r"$c_{\rm i}$")
        ax11.set_ylim([0.8 * ca_r, 1.5 * ca_t])
        ax11.set_yticks([ca_r, ca_t])
        ax11.set_yticklabels(["$c_R$", "$c_T$"])
        ax11.set_xticklabels([])
        ax11.set_xlim([0, spike_times[max_spikes - 1] + 20])
        ax11.axhline(1, ls=":", lw=1, c="C7")

        ax12.set_xlabel("$t$ / s")
        ax12.set_xlim([0, spike_times[max_spikes - 1] + 20])
        ax12.set_ylabel(r"$j_{\rm puff}$")

        ax2.set_ylim([0, 11])
        ax2.set_xlabel(r"$c_{\rm i}$")
        ax2.set_ylabel(r"$p_0(c_{\rm i})$")

        if k == 0:
            ax12.set_ylim([0, 0.7])
        else:
            ax12.set_ylim([0, 1.4])

        if k==0:
            ax2.set_xticks([ca_r, ca_t])
            ax2.set_xticklabels(["$c_R$", "$c_T$"])
            axin2.set_xticks([ca_r, ca_t])
            axin2.set_xticklabels(["$c_R$", "$c_T$"])
        else:
            ax2.set_xticks([ca_r, ca_fix, ca_t])
            ax2.set_xticklabels(["$c_R$", r"$c_{\rm i}^*$", "$c_T$"])
            axin2.set_xticks([ca_r, ca_fix, ca_t])
            axin2.set_xticklabels(["$c_R$", r"$c_{\rm i}^*$", "$c_T$"])

        if k==0:
            ax3.set_xlabel("$T$ / s")
            ax3.set_ylabel(r"$p_{\rm ISI}(T)$")

            ax3a.set_xlim([0, 75])
            ax3b.set_xlim([0, 75])

            ax3a.set_ylim([0, 0.1])
            ax3b.set_ylim([0, 0.1])
            ax3b.set_yticklabels([])
        else:
            ax6.set_xlabel("$T$ / s")
            ax6.set_ylabel(r"$p_{\rm ISI}(T)$")

            ax3a.set_xlim([0, 100])
            ax3b.set_xlim([0, 100])
            ax3a.set_ylim([0, 0.05])
            ax3b.set_ylim([0, 0.05])
            ax3a.set_xticklabels([])


    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/SUB2/figures/fig5.pdf",transparent=True)
    plt.show()
