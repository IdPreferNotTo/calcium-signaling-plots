import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    N = 10
    m = 5
    n = 4
    ip3 = 1

    xs = np.linspace(0.21, 1.1, 90)
    means_0 = []
    means_1 = []
    Ds_0 = []
    Ds_1 = []

    for x in xs:
        mean_single_0 = N*fc.mean_puff_single(x, n, m, 0.1)
        mean_single_1 = N*fc.mean_puff_single(x, n, m, 1)
        D_single_0 = np.sqrt(N)*fc.intensity_puff_single(x, n, m, 0.1)
        D_single_1 = np.sqrt(N)*fc.intensity_puff_single(x, n, m, 1)
        means_0.append(mean_single_0)
        means_1.append(mean_single_1)
        Ds_0.append(D_single_0)
        Ds_1.append(D_single_1)
    dDs_0 = []
    dDs_1 = []
    for D1, D2, in zip(Ds_0[:-1], Ds_0[1:]):
        dD = D2 - D1
        dDs_0.append(dD/0.01)
    for D1, D2, in zip(Ds_1[:-1], Ds_1[1:]):
        dD = D2 - D1
        dDs_1.append(dD/0.01)


    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    grids = gridspec.GridSpec(1, 1)
    ax0 = fig.add_subplot(grids[0])
    # ax1 = fig.add_subplot(grids[1])
    axis = [ax0]
    st.remove_top_right_axis(axis)

    ax0.grid(lw=0.5, ls="--")
    ax0.plot(xs, means_0, color=st.colors[0])
    ax0.plot(xs[:-1], dDs_0, ls=":", color=st.colors[0])
    ax0.plot(xs, means_1, label="$\mu_J(c_i)$", color=st.colors[1])
    ax0.plot(xs[:-1], dDs_1, label="$dD(c_i)/dc_i$", ls=":", color=st.colors[1])
    ax0.set_yscale("log")
    ax0.set_xlabel("$c_i$")
    ax0.legend(fancybox=False, loc=4, framealpha=1)
    ax0.axvspan(0.20, 0.33, color="C7", zorder=1, alpha=0.5)
    ax0.axvspan(1.00, 1.13, color="C7", zorder=1, alpha=0.5)
    ax0.set_xlim([0.2, 1.13])
    home = os.path.expanduser("~")
    plt.savefig(home + f"/Plots/interpretation_puff_model_N{N:d}.pdf")
    plt.show()