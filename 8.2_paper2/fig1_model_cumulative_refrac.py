import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    # Set up plot style
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3.25*1.25, 4*1.25))
    gs = gridspec.GridSpec(nrows=3, ncols=2)
    gs1 = gridspec.GridSpecFromSubplotSpec(nrows=4, ncols=1, subplot_spec=gs[0:2, 0], hspace=0.1)
    gs2 = gridspec.GridSpecFromSubplotSpec(nrows=4, ncols=1, subplot_spec=gs[0:2, 1], hspace=0.1)

    ax_ca1 = fig.add_subplot(gs1[0:2])
    ax_er1 = fig.add_subplot(gs1[2])
    ax_jp1 = fig.add_subplot(gs1[3])
    ax_ca2 = fig.add_subplot(gs2[0:2])
    ax_er2 = fig.add_subplot(gs2[2])
    ax_jp2 = fig.add_subplot(gs2[3])

    ax_isi1 = fig.add_subplot(gs[2, 0])
    ax_isi2 = fig.add_subplot(gs[2, 1])
    st.remove_top_right_axis([ax_ca1, ax_ca2, ax_er1, ax_er2, ax_jp1, ax_jp2, ax_isi1, ax_isi2])
    ax_ca1.text(0.10, 0.95, r"A$_1$", fontsize=11, transform=ax_ca1.transAxes, va='top')
    ax_ca2.text(0.10, 0.95, r"B$_1$", fontsize=11, transform=ax_ca2.transAxes, va='top')
    ax_isi1.text(0.10, 0.95, r"A$_2$", fontsize=11, transform=ax_isi1.transAxes, va='top')
    ax_isi2.text(0.10, 0.95, r"B$_2$", fontsize=11, transform=ax_isi2.transAxes, va='top')

    axis1 = [ax_ca1, ax_er1, ax_jp1, ax_isi1]
    axis2 = [ax_ca2, ax_er2, ax_jp2, ax_isi2]
    axiss = [axis1, axis2]

    ca_r = 0.2
    ca_t = 0.5
    ax_ca1.set_ylabel(r"$c_{\rm i}$")
    ax_ca1.set_yticks([ca_r, ca_t])
    ax_ca1.set_yticklabels(["$c_R$", "$c_T$"])
    ax_ca1.set_ylim([0.15, 1.5])
    ax_ca1.set_xticklabels([])

    ax_ca2.set_ylabel(r"$c_{\rm i}$")
    ax_ca2.set_yticks([0.2, 0.5])
    ax_ca2.set_yticklabels(["$c_R$", "$c_T$"])
    ax_ca2.set_ylim([0.15, 1.5])
    ax_ca2.set_xticklabels([])

    ax_er1.set_ylabel(r"$c_{\rm er}$")
    ax_er1.set_yticks([0.0, 1.0])
    ax_er1.set_xticklabels([])

    ax_er2.set_ylabel(r"$c_{\rm er}$")
    ax_er2.set_yticks([0.0, 1.0])
    ax_er2.set_xticklabels([])

    ax_jp1.set_xlabel("$t$ / s")
    ax_jp1.set_ylabel(r"$j_{\rm puff}$")

    ax_jp2.set_xlabel("$t$ / s")
    ax_jp2.set_ylabel(r"$j_{\rm puff}$")

    ax_isi1.set_xlabel("$i$")
    ax_isi1.set_ylabel("$T_i$ / s")

    ax_isi2.set_xlabel("$i$")
    ax_isi2.set_ylabel("$T_i$ / s")

    n = 5
    m = 4
    N = 10
    ip3 = 1.0

    taus = [5.0, 1.0]
    ps = [0.015, 0.060]
    eps_ers = [0.0251, 0.0794]
    tau_ers = [631.0, 63.1]
    home = os.path.expanduser("~")

    for tau, p, eps_er, tau_er, axis in zip(taus, ps, eps_ers, tau_ers, axiss):
        ax_ca = axis[0]
        ax_er = axis[1]
        ax_jp = axis[2]
        ax_isi = axis[3]

        data_traces = df.load_traces_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
        data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
        ts, cas, jpuffs, adaps = np.transpose(data_traces)
        ts_plot = []
        cas_plot = []
        as_plot = []
        jp_plot= []
        n_spikes = 0
        spike_times = []
        for t, c, jp, a in zip(ts, cas, jpuffs, adaps):
            if t >= 1_000:
                break
            ts_plot.append(t)
            as_plot.append(a)
            cas_plot.append(c)
            jp_plot.append(jp)
            if c == 0.5 and jp == 0:
                spike_times.append(t)
                ts_plot.append(t)
                cas_plot.append(0.5 + a/2)
                as_plot.append(a)
                jp_plot.append(jp)

                ts_plot.append(t)
                cas_plot.append(0.2)
                as_plot.append(a*(1 - eps_er))
                jp_plot.append(jp)
                n_spikes += 1


        for i in range(5):
            x_left = spike_times[i]
            x_right = spike_times[i+1]
            dx = spike_times[i+1] - spike_times[i]
            ax_ca.text(x_left + dx/2, 0.65, f"$T_{i+1}$", ha="center", va="center", clip_on=False)

            ax_ca.arrow(x_left + 0.05*dx, 0.55, 0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)
            ax_ca.arrow(x_right -0.05*dx, 0.55, -0.9*dx, 0, fc = "k", length_includes_head=True, head_width=0.05, head_length=5.0, lw=0.5,
                    clip_on=False)

        r0 = 1/np.mean(data_isi)
        ax_ca.plot(ts_plot, cas_plot, lw=1, color=st.colors[0])
        ax_er.plot(ts_plot, as_plot, lw=1, color=st.colors[0])
        ax_er.axhline(1/(1+eps_er*tau_er*r0))
        ax_jp.plot(ts_plot, jp_plot, lw=1, color=st.colors[0])
        ax_ca.set_xlim([0, 500])
        ax_er.set_xlim([0, 500])
        ax_jp.set_xlim([0, 500])

        idx_max = 25
        idxs = np.arange(idx_max)
        Tidxs = data_isi[:idx_max]
        popt, pcov = curve_fit(fc.exponential_Ti, idxs, Tidxs, p0=(100, 150, 2))
        T0 = popt[0]
        T8 = popt[1]
        nTr = popt[2]
        dT = T8 - T0
        Tfitixs = T0 * np.exp(-idxs / nTr) + T8 * (1 - np.exp(-idxs / nTr))
        ax_isi.set_ylim([0, 300])
        ax_isi.scatter(idxs, Tidxs, fc="w", ec=st.colors[0], s=20, zorder=3)
        ax_isi.plot(idxs, Tfitixs, lw=1, c="k", zorder=2)
        ax_isi.axvline(popt[2], lw=1, color="C7")
        ax_isi.axvspan(0, popt[2], alpha=0.5, color="C7")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig1.pdf")
    plt.show()