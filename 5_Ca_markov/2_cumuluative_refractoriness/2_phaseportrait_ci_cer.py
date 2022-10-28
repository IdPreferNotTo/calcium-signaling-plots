import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(8, 3))
    gs = gridspec.GridSpec(1, 2)
    gs1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs[0], width_ratios = [3, 1], wspace=0.0)
    gs2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs[1], width_ratios = [3, 1], wspace=0.0)

    ax1a = fig.add_subplot(gs1[0])
    ax1b = fig.add_subplot(gs1[1])
    ax2a = fig.add_subplot(gs2[0])
    ax2b = fig.add_subplot(gs2[1])
    axis = [ax1a, ax1b, ax2a, ax2b]
    st.remove_top_right_axis(axis)

    ci_min = 0.2
    ci_max = 0.5
    cer_min = 0.75
    cer_max = 1.00
    taus = [5, 1]
    ps = [0.015, 0.06]
    tau_er = 200
    eps_er = 0.1
    K = 10
    N = 5
    M = 3
    ip3 = 1.0
    cR = 0.2
    cT = 0.5
    cmap_cividis = plt.get_cmap("YlGnBu", 10)
    for tau, p, ax, axp in zip(taus, ps, [ax1a, ax2a], [ax1b, ax2b]):
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
        ax.hist2d(cis, cers, bins=[25,25], alpha=0.5, cmap=cmap_cividis, zorder=0)
        data_exit = [data for data in data_ci if data[1] == 0.5 and data[2] == 0]
        data_exit = np.asarray(data_exit)
        axp.hist(data_exit[:,3], density=True, color=st.colors[0], bins=25, orientation="horizontal")

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
        ax.set_ylim([cer_min, cer_max])
        axp.set_xlabel(r"$p(c_{\rm er}; c_T)$")
        axp.set_ylim([cer_min, cer_max])
        axp.set_yticks([])


    plt.show()