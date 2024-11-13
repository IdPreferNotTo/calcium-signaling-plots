import os

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 1.35 * 1.25
    fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(w, h))
    ax1 = axs[0]
    ax2 = axs[1]


    st.remove_top_right_axis([ax1, ax2])
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')

    ax1.set_ylabel(r"$T_i$ / s")
    ax1.set_xlabel(r"$i$")
    ax2.set_ylabel(r"$T_i$ / s")
    ax2.set_xlabel(r"$i$")

    #------------------ EXAMPLES ---------------------------------------------------------------------------------------
    home = os.path.expanduser("~")
    ax1.set_xlim([0, 20])
    ax1.set_xticks([0, 5, 10, 15, 20])
    ax1.set_ylim([0, 400])
    ax2.set_xlim([0, 20])
    ax2.set_xticks([0, 5, 10, 15, 20])
    ax2.set_ylim([0, 400])

    pars = [[4.42, 1.72e-02, 4.83e+02, 8.02e-02], [9.92, 1.04e-02, 1.04e+03, 3.13e-02]]
    for i, (ax, par) in enumerate(zip([ax1, ax2], pars)):
        data_isi = df.load_spike_times_markov_transient(par[0], par[1], par[2], par[3])
        rows, cols = data_isi.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs]  # calculate the mean column wise
        popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        t0 = means_Tidx[0]
        t8 = means_Tidx[-1]
        func = lambda x, ntr: fc.exponential_Ti(x, t0, t8, ntr)  # fix t0, t8
        popt, pcov = curve_fit(func, idxs, means_Tidx, p0=[2.])
        ntr = popt[0]
        Tis = [t0 * np.exp(-idx / ntr) + t8 * (1 - np.exp(-idx / ntr)) for idx in idxs]

        ax.plot(idxs, Tis, c="k", lw=1, label=rf"$n_{{\rm tr}} = {ntr:.1f}$" + "\n"
                                                 + rf"$\Delta T = {round(t8 - t0,-1):.0f}$s")
        ax.scatter(idxs, means_Tidx, s=15, fc="w", ec=st.colors[1])
        ax.axhline(t0, ls=":", c="C7")
        ax.axhline(t8, ls=":", c="C7")
        ax.axvspan(0, math.ceil(2*ntr), facecolor="C7", alpha=0.3, zorder=0)

        if i == 0:
            leg = ax.legend(fancybox=False, framealpha=1., loc=4, edgecolor="k", fontsize=7)
            leg.get_frame().set_linewidth(0.75)
            for t in leg.get_texts():
                t.set_ha("center")
        else:
            leg = ax.legend(fancybox=False, framealpha=1., loc=1, edgecolor="k", fontsize=7)
            leg.get_frame().set_linewidth(0.75)
            for t in leg.get_texts():
                t.set_ha("center")

    plt.savefig(home + f"/Desktop/model_fit.png", dpi=300)
    plt.show()