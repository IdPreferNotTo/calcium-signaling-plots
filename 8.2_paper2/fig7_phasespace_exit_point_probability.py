import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 3))
    gs = gridspec.GridSpec(2, 3)
    gs1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs[0], width_ratios = [3, 1], wspace=0.0)
    gs2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs[1], width_ratios = [3, 1], wspace=0.0)
    gs3 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs[2], width_ratios=[3, 1], wspace=0.0)

    ax1a = fig.add_subplot(gs1[0])
    ax1b = fig.add_subplot(gs1[1])
    ax1c = fig.add_subplot(gs[1, 0])
    ax2a = fig.add_subplot(gs2[0])
    ax2b = fig.add_subplot(gs2[1])
    ax2c = fig.add_subplot(gs[1, 1])
    ax3a = fig.add_subplot(gs3[0])
    ax3b = fig.add_subplot(gs3[1])
    ax3c = fig.add_subplot(gs[1, 2])

    axis = [ax1a, ax1b, ax1c, ax2a, ax2b, ax2c, ax3a, ax3b, ax3c]
    st.remove_top_right_axis(axis)

    ci_min = 0.2
    ci_max = 0.5
    cer_min = 0.50
    cer_max = 1.00
    taus = [5, 5, 1]
    ps = [0.015, 0.015, 0.06]
    tau_ers = [501, 501, 501]
    eps_ers = [0.0398, 0.398, 0.0398]
    K = 10
    N = 5
    M = 3
    ip3 = 1.0
    cR = 0.2
    cT = 0.5
    cmap_cividis = plt.get_cmap("YlGnBu", 10)
    for tau, p, tau_er, eps_er, ax, axp, axisi in zip(taus, ps, tau_ers, eps_ers, [ax1a, ax2a, ax3a], [ax1b, ax2b, ax3b], [ax1c, ax2c, ax3c]):
        ax.axvline(cR, ls=":", lw=1, c="k")
        ax.axvline(cT, ls="--", lw=1, c="k")

        cis = np.linspace(ci_min, ci_max, 51, endpoint=True)
        cers = np.linspace(cer_min, cer_max, 51, endpoint=True)
        dxs = np.zeros((51, 51))
        dys = np.zeros((51, 51))
        dvs = np.zeros((51, 51))
        for i, cer in enumerate(cers):
            for j, ci in enumerate(cis):
                if ci == 0:
                    dxs[i, j] = -(ci - cR) / tau
                else:
                    dxs[i,j] = -(ci - cR) / tau + cer * p * K * fc.mean_jp_single_theory(ci, N, M, ip3)
                dys[i,j] = (1 - cer) / tau_er
                dvs[i,j] = np.sqrt(dxs[i, j]**2 + dys[i, j]**2)
        strm = ax.streamplot(cis, cers, dxs, dys, color="C7", arrowsize=0.75, linewidth=0.3)
        data_ci = df.load_traces_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
        ts, cis, jps, cers = np.transpose(data_ci)
        pbif = fc.pbif(tau, K=10, N=5, M=3, s=1.0)
        cer_crit = pbif/p
        data_exit = [data for data in data_ci if data[1] == 0.5 and data[2] == 0]
        data_exit_excitable = [data for data in data_ci if data[1] == 0.5 and data[2] == 0 and data[3] < cer_crit]
        data_exit_meandriven = [data for data in data_ci if data[1] == 0.5 and data[2] == 0 and data[3] > cer_crit]
        data_exit = np.asarray(data_exit)
        numbers, bins, patches = axp.hist(data_exit[5:,3], density=True, alpha=.75, color=st.colors[0], bins=50, orientation="horizontal")
        for idx, bin in enumerate(bins):
            if bin > cer_crit:
                break
        for i in range(idx):
            patches[i].set_facecolor(st.colors[1])
        for i in range(idx, len(patches)):
            patches[i].set_facecolor(st.colors[5])

        ax.hist2d(cis, cers, bins=[25,25], alpha=0.5, cmap=cmap_cividis, zorder=0)

        data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)
        mean_isi = np.mean(data_isi[10:])
        times, cis, jps, cers = np.transpose(data_exit)
        isi_excitable = []
        isi_mean_driven = []
        for isi, cer in zip(data_isi, cers):
            if cer > cer_crit:
                isi_mean_driven.append(isi)
            else:
                isi_excitable.append(isi)

        axisi.hist([isi_mean_driven, isi_excitable], bins=25, color=[st.colors[5], st.colors[1]], stacked=True,
                 alpha=0.75)

        if tau == 5:
            x_lc = []
            y_lc = []
            ci = 0.2
            cer = 1.00
            dt = 0.01
            spike_count = 0
            while spike_count < 100:
                ci += -(ci - cR) / tau + cer * p * K * fc.mean_jp_single_theory(ci, N, M, ip3)
                cer += (1 - cer) / tau_er
                if ci > 0.5:
                    ci = cR
                    cer -= eps_er * cer
                    spike_count += 1
            spike_count = 0
            while spike_count < 1:
                x_lc.append(ci)
                y_lc.append(cer)
                ci += -(ci - cR) / tau + cer * p * K * fc.mean_jp_single_theory(ci, N, M, ip3)
                cer += (1 - cer) / tau_er
                if ci > 0.5:
                    x_lc.append(cT)
                    y_lc.append(cer)
                    cer -= eps_er * cer
                    x_lc.append(cT)
                    y_lc.append(cer)
                    ci = cR
                    x_lc.append(ci)
                    y_lc.append(cer)
                    spike_count += 1
            ax.plot(x_lc, y_lc, lw=1, c="k")

        ax.set_xlabel(r"$c_{\rm i}$")
        ax.set_ylabel(r"$c_{\rm er}$")
        ax.set_xlim([ci_min, ci_max])
        ax.set_xticks([cR, cT])
        ax.set_xticklabels(["$c_R$", "$c_T$"])
        axp.set_xlabel(r"$p_{\rm exit}(c_{\rm er})$")
        axp.set_yticks([])
        axp.set_xticks([])
        axisi.set_xlabel(r"$t$")
        axisi.set_ylabel(r"$p_{\rm ISI}(t)$")

        cerminus = 1 - eps_er / (1 - (1 - eps_er)*np.exp(-mean_isi/tau_er))
        cermean = 1 / (1 + eps_er*tau_er/mean_isi)
        Tdet = tau_er*np.log((1-cerminus)/(1-cermean))
        print(Tdet)
        axp.axhline(cerminus, c="k")
        axp.axhline(cermean, c="k")
        axisi.axvline(Tdet, c="C7", ls=":")
        data_isi_no_adap = df.load_spike_times_markov(tau, p, cer=False)
        mean_isi_no_adap = np.mean(data_isi_no_adap)
        axisi.hist(data_isi_no_adap, alpha=0.5, bins=25, color="k")

    ax1a.set_ylim([0.8, cer_max])
    ax1b.set_ylim([0.8, cer_max])
    ax2a.set_ylim([0.5, cer_max])
    ax2b.set_ylim([0.5, cer_max])
    ax3a.set_ylim([0.8, cer_max])
    ax3b.set_ylim([0.8, cer_max])
    home = os.path.expanduser("~")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig7.pdf", transparent=True)
    plt.show()