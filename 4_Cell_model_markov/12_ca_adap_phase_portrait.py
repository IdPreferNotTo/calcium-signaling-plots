import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import functions as fc
import styles as st
if __name__ == "__main__":
    tau = 0.2
    jca = 0.501
    taua = 500
    ampa = 0.05
    N = 10
    n = 5
    m= 4
    ip3 = 1.0

    ca_r = 0.33
    ca_t = 1.00
    ca_i_min = 0.2
    ca_er_min = 0.5
    ca_is = np.linspace(ca_i_min, 1, 51, endpoint=True)
    ca_ers = np.linspace(ca_er_min, 1, 51, endpoint=True)
    dxs = np.zeros((51, 51))
    dys = np.zeros((51, 51))
    dvs = np.zeros((51, 51))
    for i, ca_er in enumerate(ca_ers):
        for j, ca_i in enumerate(ca_is):
            if ca_i == 0:
                dxs[i, j] = -(ca_i - ca_r) / tau
            else:
                dxs[i,j] = -(ca_i - ca_r) / tau + ca_er * jca * N * fc.mean_puff_single(ca_i, n, m, ip3)
            dys[i,j] = (1 - ca_er) / taua
            dvs[i,j] = np.sqrt(dxs[i, j]**2 + dys[i, j]**2)

    ca_i_fix = tau*jca * N * fc.mean_puff_single(ca_i, n, m, ip3) + ca_r
    x_lc = []
    y_lc = []
    ca_i = 0.33
    ca_er = 1.00
    dt = 0.01
    spike_count = 0
    while spike_count < 100:
        ca_i += -(ca_i - 0.33) / tau + ca_er * jca * N * fc.mean_puff_single(ca_i, n, m, ip3)
        ca_er += (1 - ca_er) / taua
        if ca_i > 1.:
            ca_i = 0.33
            ca_er -= ampa * ca_er
            spike_count += 1
    spike_count = 0
    while spike_count < 1:
        x_lc.append(ca_i)
        y_lc.append(ca_er)
        ca_i += -(ca_i - 0.33) / tau + ca_er * jca * N * fc.mean_puff_single(ca_i, n, m, ip3)
        ca_er += (1 - ca_er) / taua
        if ca_i > 1.:
            x_lc.append(1)
            y_lc.append(ca_er)
            ca_er -= 0.2 * ca_er
            x_lc.append(1)
            y_lc.append(ca_er)
            ca_i = 0.33
            x_lc.append(ca_i)
            y_lc.append(ca_er)
            spike_count += 1

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    ax.set_ylabel(r"[Ca\textsuperscript{2+}]$_{\rm ER}$ [a.u.]")
    ax.set_xlabel(r"[Ca\textsuperscript{2+}]$_{\rm i}$ [a.u.]")
    ax.set_xlim([xmin, 1.05])
    ax.set_ylim([ymin, 1.05])
    st.remove_top_right_axis([ax])
    #strm = ax.streamplot(xs, ys, dxs, dys, color=dvs, arrowsize=0.75, linewidth=0.3, cmap="cividis")
    strm = ax.streamplot(ca_is, ca_ers, dxs, dys, color="C7", arrowsize=0.75, linewidth=0.3)

    ax.plot(x_lc, y_lc, lw=1, c="k")
    ax.arrow(1, 0.75, dx=0, dy=-0.05, shape="full", fc="k", lw=0, length_includes_head=True,
              head_length=0.05, head_width=.025)
    ax.arrow((1 + 0.33)/2, y_lc[-1], dx=-0.05, dy=0.0, shape="full", fc="k", lw=0, length_includes_head=True,
              head_length=.05, head_width=.025)

    ax.annotate("depletion", xy=(1 + 0.05, (y_lc[-3] + y_lc[-2])/2), va="center", rotation=-90)
    ax.annotate("reset", xy=(0.6, y_lc[-1] - 0.05), va="center", backgroundcolor="white")

    ax.axvline(0.33, ls=":", lw=1, c="k")
    ax.axvline(1, ls="--", lw=1, c="k")
    home = os.path.expanduser("~")
    plt.savefig(home + f"/Plots/8_markov_ca_adap_phaseportrait_tau{tau:.2e}j{jca:.2e}.pdf", transparent=True)
    plt.show()