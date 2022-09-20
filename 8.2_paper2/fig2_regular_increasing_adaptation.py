import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st

def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))

def inverse_gaussian(T, CV):
    ps = []
    ts = np.linspace(0, 2*T, 500)
    for t in ts[1:]:
        p = np.sqrt(T / (2 * np.pi * (CV**2) * (t ** 3))) * np.exp(-(t - T) ** 2 / (2 * T * (CV**2) * t))
        ps.append(p)
    return ts[1:], ps

if __name__ == "__main__":
    ampa1 = 0
    ampa2 = 0.1 #np.logspace(-2, 1, 50)[10]
    ampa3 = 0.5 #np.logspace(-2, 1, 50)[15]
    taua = 100 #829
    tau = 10. # 10.5
    j = 0.015 #0.0146
    amps = [ampa1, ampa2, ampa3]

    home = os.path.expanduser("~")

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(9, 6))
    gs = gridspec.GridSpec(nrows=4, ncols=3)
    gs_a1 = gridspec.GridSpecFromSubplotSpec(nrows=4, ncols=1, subplot_spec=gs[0:2, 0], hspace=0.05)
    gs_b1 = gridspec.GridSpecFromSubplotSpec(nrows=4, ncols=1, subplot_spec=gs[0:2, 1], hspace=0.05)
    gs_c1 = gridspec.GridSpecFromSubplotSpec(nrows=4, ncols=1, subplot_spec=gs[0:2, 2], hspace=0.05)

    ax_a11 = fig.add_subplot(gs_a1[0:2])
    ax_a12 = fig.add_subplot(gs_a1[2])
    ax_a13 = fig.add_subplot(gs_a1[3])
    ax_a2 = fig.add_subplot(gs[2, 0])
    ax_a3 = fig.add_subplot(gs[3, 0])
    axis_a = [ax_a11, ax_a12, ax_a13, ax_a2, ax_a3]

    ax_b11 = fig.add_subplot(gs_b1[0:2])
    ax_b12 = fig.add_subplot(gs_b1[2])
    ax_b13 = fig.add_subplot(gs_b1[3])
    ax_b2 = fig.add_subplot(gs[2, 1])
    ax_b3 = fig.add_subplot(gs[3, 1])
    axis_b = [ax_b11, ax_b12, ax_b13, ax_b2, ax_b3]

    ax_c11 = fig.add_subplot(gs_c1[0:2])
    ax_c12 = fig.add_subplot(gs_c1[2])
    ax_c13 = fig.add_subplot(gs_c1[3])
    ax_c2 = fig.add_subplot(gs[2, 2])
    ax_c3 = fig.add_subplot(gs[3, 2])
    axis_c = [ax_c11, ax_c12, ax_c13, ax_c2, ax_c3]

    st.remove_top_right_axis(axis_a + axis_b + axis_c)

    ax_a11.text(0.05, 0.95, "A$_{i}$", fontsize=13, transform=ax_a11.transAxes, va='top')
    ax_a2.text(0.05, 0.95, "A$_{ii}$", fontsize=13, transform=ax_a2.transAxes, va='top')
    ax_a3.text(0.05, 0.95, "A$_{iii}$", fontsize=13, transform=ax_a3.transAxes, va='top')

    ax_b11.text(0.05, 0.95, "B$_{i}$", fontsize=13, transform=ax_b11.transAxes, va='top')
    ax_b2.text(0.05, 0.95, "B$_{ii}$", fontsize=13, transform=ax_b2.transAxes, va='top')
    ax_b3.text(0.05, 0.95, "B$_{iii}$", fontsize=13, transform=ax_b3.transAxes, va='top')

    ax_c11.text(0.05, 0.95, "C$_{i}$", fontsize=13, transform=ax_c11.transAxes, va='top')
    ax_c2.text(0.05, 0.95, "C$_{ii}$", fontsize=13, transform=ax_c2.transAxes, va='top')
    ax_c3.text(0.05, 0.95, "C$_{iii}$", fontsize=13, transform=ax_c3.transAxes, va='top')

    for axis, amp in zip([axis_a, axis_b, axis_c], amps):
        ax_11, ax_12, ax_13, ax_2, ax_3 = axis
        ax_11.text(0.5, 0.95, f"$\epsilon = {amp:.2f}$", fontsize=13, transform=ax_11.transAxes, ha='center', va='top')
        folder = home + "/Data/calcium_spikes_markov/Data_adap"
        file_ca = f"/ca_markov_ip1.00_taua{taua:.2e}_ampa{amp:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
        file_isi = f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{amp:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
        data_ca = np.loadtxt(folder + file_ca)
        ts, cas, jpuffs, adaps = np.transpose(data_ca)
        isis = np.loadtxt(folder + file_isi)

        count = 0
        spike_times = []
        max_i = 0
        max_spikes = 0
        sum_isi = 0
        while sum_isi < 300:
            sum_isi += isis[max_spikes]
            max_spikes +=1

        while count < max_spikes:
            if cas[max_i] == 1:
                count += 1
                spike_times.append(ts[max_i])
            max_i += 1

        ax_11.plot(ts[:max_i], cas[:max_i])
        ax_12.plot(ts[:max_i], adaps[:max_i])
        ax_12.set_ylabel(r"$\Delta c_{er}$")
        ax_12.set_ylim([0.4, 1.1])
        ax_12.set_yticks([0.5, 1.0])
        ax_12.set_xticklabels([])
        ax_12.set_xlim([0, 300])

        ax_13.plot(ts[:max_i], jpuffs[:max_i])
        ax_13.set_ylabel(r"$j_{\rm puff}$")
        ax_13.set_xlabel("$t$ / s")
        ax_13.set_xlim([0, 300])

        for i in range(max_spikes-1):
            if i == 0:
                x_left = 0
                x_right = spike_times[i]
                dx = x_right
            else:
                x_left = spike_times[i-1]
                x_right = spike_times[i]
                dx = x_right - x_left
            if i < 4:
                ax_11.text(x_left + dx/2, 1.20, f"$T_{i}$", ha="center", va="center", clip_on=False)
            if i == 4:
                ax_11.text(x_left + dx / 2, 1.20, f"$...$", ha="center", va="center", clip_on=False)

            ax_11.arrow(x_left + 0.05*dx, 1.10, 0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)
            ax_11.arrow(x_right -0.05*dx, 1.10, -0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)

        #ax0.set_xlabel("$t$ / s")
        ca_t = 1.
        ca_r = 1/3
        ax_11.set_ylabel(r"$c_i$")
        ax_11.set_ylim([0.8 * ca_r, 1.8 * ca_t])
        ax_11.set_yticks([ca_r, ca_t])
        ax_11.set_yticklabels(["$c_R$", "$c_T$"])
        ax_11.set_xticklabels([])
        ax_11.set_xlim([0, 300])
        ax_11.axhline(1, ls=":", lw=1, c="C7")

        nr_isis = 10
        ax_2.set_xlim([0, nr_isis])
        #ax_2.set_yticks([0, 50, 100, 150])
        ax_2.set_ylim([0, 150])
        ax_2.set_xlabel("$i$")
        ax_2.set_ylabel("$T_i$ / s")
        ax_2.scatter(np.arange(0, nr_isis), isis[0:nr_isis], fc="w", ec=st.colors[1], s=20, zorder=3)
        if amp == 0:
            ax_2.text(nr_isis / 2, 2*np.mean(isis), "$T_0 = T_\infty$", ha="center")
            ax_2.axhline(np.mean(isis), ls=":", lw=1, c="k")
        else:
            index_isis = np.arange(len(isis))
            popt, pcov = curve_fit(transient_func, index_isis, isis, p0=(100, 150, 2))
            ax_2.axvspan(0, popt[2], alpha=0.3, color="C7")

            print(popt)
            ISI_fit = popt[0] * np.exp(-index_isis / popt[2]) + popt[1] * (1. - np.exp(-index_isis / popt[2]))
            ax_2.plot(index_isis, ISI_fit, lw=1, c="k")
            ax_2.axhline(popt[0], ls=":", lw=1, c="k")
            ax_2.axhline(popt[1], ls=":", lw=1, c="k")
            ax_2.text(nr_isis / 2, popt[0] * 1.5, "$T_0$", ha="center")
            ax_2.text(nr_isis / 2, popt[1] * 1.2, "$T_\infty$", ha="center")

        mean_isis = np.mean(isis)
        std_isis = np.std(isis)
        cv_isis = std_isis/mean_isis
        ts_inv_gaussian, ps_inv_gaussian = inverse_gaussian(mean_isis, cv_isis)
        ax_3.set_xlim(0, 2*mean_isis)
        ax_3.set_xticks([0, mean_isis, 2*mean_isis])
        ax_3.set_xticklabels(["$0$", r"$\langle T \rangle$", r"2$\langle T \rangle$"])
        ax_3.set_xlabel("$T$")
        ax_3.set_ylabel("$P(T)$")
        ax_3.hist(isis, bins=20, color=st.colors[1], density=True, alpha=0.6)
        ax_3.plot(ts_inv_gaussian, ps_inv_gaussian, color="k",  ls="--")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig2.pdf", transparent=True)
    plt.show()