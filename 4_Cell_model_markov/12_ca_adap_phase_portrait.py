import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from functions import *

if __name__ == "__main__":
    tau = 2.81
    jca = 0.0728
    taua = 100
    ampa = 0.2
    n_cl = 10

    xmin = 0.2
    ymin = 0.5
    xs = np.linspace(xmin, 1, 21, endpoint=True)
    ys = np.linspace(ymin, 1, 21, endpoint=True)
    dxs = np.zeros((21, 21))
    dys = np.zeros((21, 21))
    dvs = np.zeros((21, 21))
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            if x == 0:
                dxs[i, j] = -(x - 0.33)/tau
            else:
                dxs[i,j] = -(x - 0.33)/tau + y*jca*n_cl*mean_puff_single(x)
            dys[i,j] = (1-y)/taua
            dvs[i,j] = np.sqrt(dxs[i, j]**2 + dys[i, j]**2)
    #dvs = np.sqrt(dxss ** 2 + dyss ** 2)

    x_lc = []
    y_lc = []
    x = 0.33
    y = 1.00
    dt = 0.01
    spike_count = 0
    while spike_count < 100:
        x += -(x - 0.33)/tau + y*jca*n_cl*mean_puff_single(x)
        y += (1-y)/taua
        if x > 1.:
            x = 0.33
            y -= 0.2*y
            spike_count += 1
    spike_count = 0
    while spike_count < 1:
        x_lc.append(x)
        y_lc.append(y)
        x += -(x - 0.33) / tau + y * jca * n_cl * mean_puff_single(x)
        y += (1 - y) / taua
        if x > 1.:
            x_lc.append(1)
            y_lc.append(y)
            y -= 0.2 * y
            x_lc.append(1)
            y_lc.append(y)
            x = 0.33
            x_lc.append(x)
            y_lc.append(y)
            spike_count += 1

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    ax.set_ylabel(r"[Ca\textsuperscript{2+}]$_{\rm ER}$ [a.u.]")
    ax.set_xlabel(r"[Ca\textsuperscript{2+}]$_{\rm i}$ [a.u.]")
    ax.set_xlim([xmin, 1.05])
    ax.set_ylim([ymin, 1.05])
    remove_top_right_axis([ax])
    #strm = ax.streamplot(xs, ys, dxs, dys, color=dvs, arrowsize=0.75, linewidth=0.3, cmap="cividis")
    strm = ax.streamplot(xs, ys, dxs, dys, color="C7", arrowsize=0.75, linewidth=0.3)

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
    plt.savefig(home + f"/Data/Calcium/Plots/8_markov_ca_adap_phaseportrait_tau{tau:.2e}j{jca:.2e}.pdf", transparent=True)
    plt.show()