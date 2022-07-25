import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    color = st.Colors
    fig = plt.figure(tight_layout=True, figsize=(9, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    gs11 = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs[0, 0], hspace=0.1)
    gs12 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0, 1], wspace=0.1)
    gs13 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0, 2], wspace=0.1)

    gs21 = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs[1, 0], hspace=0.1)
    gs22 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[1, 1], wspace=0.1)
    gs23 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[1, 2], wspace=0.1)

    ax11= fig.add_subplot(gs11[0:2])
    ax12 = fig.add_subplot(gs11[2])
    ax2 = fig.add_subplot(gs12[0])
    axin2 = ax2.inset_axes([0.65, 0.65, 0.35, 0.35])
    ax3 = fig.add_subplot(gs13[0])

    ax41 = fig.add_subplot(gs21[0:2])
    ax42 = fig.add_subplot(gs21[2])
    ax5 = fig.add_subplot(gs22[0])
    axin5 = ax5.inset_axes([0.22, 0.65, 0.35, 0.35])
    ax6 = fig.add_subplot(gs23[0])


    st.remove_top_right_axis([ax11, ax12, ax2, axin2, ax3, ax41, ax42, ax5, axin5, ax6])
    axis1 = [ax11, ax12, ax2, axin2, ax3]
    axis2 = [ax41, ax42, ax5, axin5, ax6]
    axiss = [axis1, axis2]

    ax11.text(0.05, 0.95, r"A$_{\rm i}$", fontsize=13, transform=ax11.transAxes, va='top')
    ax2.text(0.05, 0.95, r"A$_{\rm ii}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"A$_{\rm iii}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax41.text(0.05, 0.95, r"B$_{\rm i}$", fontsize=13, transform=ax41.transAxes, va='top')
    ax5.text(0.05, 0.95, r"B$_{\rm ii}$", fontsize=13, transform=ax5.transAxes, va='top')
    ax6.text(0.05, 0.95, r"B$_{\rm iii}$", fontsize=13, transform=ax6.transAxes, va='top')
    home = os.path.expanduser("~")
    # Plot regular spiketrain
    taus = [4.47, 1.26]
    js = [0.0178, 0.0562]
    for k, (tau, j, axis) in enumerate(zip(taus, js, axiss)):

        ax11, ax12, ax2, axin2, ax3 = axis
        ca_r = 0.33
        ca_t = 1.00
        N = 10
        n = 5
        m = 3
        ip3 = 1.0
        r_opn_single = 0.2
        r_ref = 5

        # Get Data
        folder = home + "/Data/calcium_spikes_markov/Data_no_adap/"
        file_calcium = f"ca_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        data_calcium = np.loadtxt(folder + file_calcium)
        data_isis = np.loadtxt(folder + file_spikes)
        ts, cas, jpuffs, adaps = np.transpose(data_calcium)

        folder_langevin = home + "/Data/calcium_spikes_langevin_strat/Data_no_adap/"
        file_spikes_langevin = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_0.dat"
        data_isis_langevin = np.loadtxt(folder_langevin + file_spikes_langevin)
        mean_isi_langevin_sim = np.mean(data_isis_langevin)
        std_isi_langevin_sim = np.std(data_isis_langevin)
        cv_isi_langevin_sim = std_isi_langevin_sim/mean_isi_langevin_sim

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
        cv = np.loadtxt(file_cv_theory)
        cv_isi_langevin = cv
        if k == 0:
            plt_color = color.palette[1]
        else:
            plt_color = color.orange[1]
            
        #Get drift
        xs = np.linspace(0.3, 1, 100)
        dx = xs[1] - xs[0]
        ca_drifts = []
        for x in xs:
            ca_drift = -(x- ca_r)/tau + j*N*fc.mean_puff_single(x, n, m, ip3, r_opn_single, r_ref)
            ca_drifts.append(20*ca_drift)

        # Plots ax11 ax12
        count = 0
        spike_times = []
        max_i = 0
        max_spikes = 5
        while count < max_spikes:
            if cas[max_i] == 1:
                count += 1
                spike_times.append(ts[max_i])
            max_i += 1
        ax11.plot(ts, cas, color=plt_color)
        ax12.plot(ts_puff, [j*n for n in ns_puff], color=plt_color)
        #ax12.plot(ts[:max_i], jpuffs[:max_i], color=plt_color)

        #Plot ax2
        cas = [ca for ca in cas if ca !=1]
        ax2.hist(cas, bins=25, density=True, alpha=0.5, color=plt_color)
        ax2.plot(ca_theory[::100], p0_theory[::100], color=color.palette[5])
        axin2.plot(xs, ca_drifts, color=plt_color)
        axin2.set_title(r"$g(c_{\rm i})$")
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
            ax3.hist(data_isis, bins=20, color=color.palette[1], density=True, alpha=0.5)
        else:
            ax3.hist(data_isis, bins=50, color=color.orange[1], density=True, alpha=0.5)


        ts_inv_gau = np.linspace(1, 200, 1001)
        ts_exp = np.linspace(1, 200, 1001)
        inv_gaus_dis = []
        exp_dis = []
        for t in ts_inv_gau:
            p = np.sqrt(mean_isi / (2 * np.pi * np.power(cv_isi, 2) * (t ** 3))) * np.exp(
                -(t - mean_isi) ** 2 / (2 * mean_isi * np.power(cv_isi, 2) * t))
            inv_gaus_dis.append(p)
        ax3.plot(ts_inv_gau, inv_gaus_dis, color=color.palette[5], ls="--", label="Two-component")

        ts_inv_gau = np.linspace(1, 200, 1001)
        ts_exp = np.linspace(1, 200, 1001)
        inv_gaus_dis = []
        exp_dis = []
        for t in ts_inv_gau:
            p = np.sqrt(mean_isi_langevin / (2 * np.pi * np.power(cv_isi_langevin_sim, 2) * (t ** 3))) * np.exp(
                -(t - mean_isi_langevin) ** 2 / (2 * mean_isi_langevin * np.power(cv_isi_langevin_sim, 2) * t))
            inv_gaus_dis.append(p)
        ax3.plot(ts_inv_gau, inv_gaus_dis, color=color.palette[5], ls=":", label="Langevin")
        if k == 0:
            ax3.legend(fancybox=False, fontsize=9, loc=1, title="Inverse Gaussian")

        for i in range(max_spikes-1):
            x_left = spike_times[i]
            x_right = spike_times[i+1]
            dx = spike_times[i+1] - spike_times[i]
            ax11.text(x_left + dx/2, 1.20, f"$T_{i+1}$", ha="center", va="center", clip_on=False)

            ax11.arrow(x_left + 0.05*dx, 1.10, 0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)
            ax11.arrow(x_right -0.05*dx, 1.10, -0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)

        #ax0.set_xlabel("$t$ / s")
        ax11.set_ylabel(r"$c_{\rm i}$")
        ax11.set_ylim([0.8 * ca_r, 1.5 * ca_t])
        ax11.set_yticks([ca_r, ca_t])
        ax11.set_yticklabels(["$c_0$", "$c_T$"])
        ax11.set_xticklabels([])
        ax11.set_xlim([0, spike_times[max_spikes - 1] + 20])
        ax11.axhline(1, ls=":", lw=1, c="C7")

        ax12.set_xlabel("$t$ / s")
        ax12.set_xlim([0, spike_times[max_spikes - 1] + 20])
        ax12.set_ylabel(r"$j_{\rm puff}$")
        ax12.set_ylim([0,1.5])

        ax2.set_ylim([0, 4])
        ax2.set_xlabel(r"$c_{\rm i}$")
        ax2.set_ylabel(r"$p_0(c_{\rm i})$")

        ax2.set_xticks([0.33, 1.0])
        ax2.set_xticklabels(["$c_0$", "$c_T$"])
        axin2.set_xticks([ca_r, ca_t])
        axin2.set_xticklabels(["$c_0$", "$c_T$"])


        #ax3.text(0.6, 0.90, rf"$\langle T \rangle = {mean_isi:.1f}$ s", transform=ax3.transAxes, va='top')
        #ax3.text(0.61, 0.75, rf"$CV_T = {cv_isi:.2f}$", transform=ax3.transAxes, va='top')
        ax3.set_xlabel("$T$ / s")
        ax3.set_ylabel(r"$p_{\rm ISI}(T)$")
        ax3.set_xlim([0, 200])
        ax3.set_ylim([0, 0.06])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig6.pdf",transparent=True)
    plt.show()
