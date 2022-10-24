import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 3))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    axis = [ax1, ax2]
    st.remove_top_right_axis(axis)

    ci_min = 0.0
    ci_max = 0.6
    cer_min = 0.5
    cer_max = 1.0
    taus = [5, 1]
    ps = [0.015, 0.06]
    tau_er = 100
    eps_er = 0.1
    K = 10
    N = 5
    M = 3
    ip3 = 1.0
    cR = 0.2
    cT = 0.5

    for tau, p, ax in zip(taus, ps, [ax1, ax2]):
        ax.set_xlabel(r"$c_{\rm i}$")
        ax.set_ylabel(r"$c_{\rm er}$")
        ax.set_xlim([ci_min, ci_max])
        ax.set_ylim([cer_min, cer_max])
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
            #strm = ax.streamplot(xs, ys, dxs, dys, color=dvs, arrowsize=0.75, linewidth=0.3, cmap="cividis")
            ax.plot(x_lc, y_lc, lw=1, c="k")
        #ax.arrow(1, 0.75, dx=0, dy=-0.05, shape="full", fc="k", lw=0, length_includes_head=True, head_length=0.05, head_width=.025)
        #ax.arrow((1 + 0.33)/2, y_lc[-1], dx=-0.05, dy=0.0, shape="full", fc="k", lw=0, length_includes_head=True, head_length=.05, head_width=.025)
        #ax.annotate("depletion", xy=(1 + 0.05, (y_lc[-3] + y_lc[-2])/2), va="center", rotation=-90)
        #ax.annotate("reset", xy=(0.6, y_lc[-1] - 0.05), va="center", backgroundcolor="white")


    plt.show()