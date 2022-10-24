import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axis = [ax1, ax3, ax2, ax4]
    axis_idx = [ax1, ax3]
    axis_return = [ax2, ax4]
    st.remove_top_right_axis(axis)

    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    for tau, p, ax_idx, ax_ret in zip(taus, ps, axis_idx, axis_return):
        tau_er = 500
        eps_er = 0.025
        data_isi = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er)
        rows, cols = data_isi.shape
        idx_max = cols
        idxs = np.arange(idx_max)
        means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise

        for row in np.arange(1):
            ax_idx.scatter(idxs, data_isi[row, :], fc="w", ec=st.colors[0], alpha=0.5, s=20, zorder=3)
        ax_idx.scatter(idxs, means_Tidx, fc="w", ec=st.colors[5], s=20, zorder=3)
        popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=(100, 150, 2))
        ax_idx.axvspan(0, popt[2], alpha=0.3, color="C7")

        # print(popt)
        ISI_fit = popt[0] * np.exp(-idxs / popt[2]) + popt[1] * (1. - np.exp(-idxs / popt[2]))
        ax_idx.plot(idxs, ISI_fit, lw=1, c="k")
        ax_idx.axhline(popt[0], ls=":", lw=1, c="k")
        ax_idx.axhline(popt[1], ls=":", lw=1, c="k")
        if int(popt[2]) != 0:
            ax_idx.set_xlim([0, 5*int(popt[2])])
        else:
            ax_idx.set_xlim(([0, 5]))

        for row in np.arange(1):
            ax_ret.scatter(np.asarray(popt[1])-data_isi[row, :-1], np.asarray(popt[1])-data_isi[row, 1:], fc="w", ec=st.colors[0], alpha=0.5, s=20, zorder=3)
        ax_ret.scatter(np.asarray(popt[1])-means_Tidx[:-1], np.asarray(popt[1])-means_Tidx[1:], fc="w", ec=st.colors[5], s=20, zorder=3)
        ax_ret.axhline(0, ls=":", lw=1, c="k")
        ax_ret.axvline(0, ls=":", lw=1, c="k")
        xs = np.linspace(0, popt[1])
        ys = np.exp(-1/popt[2]) * xs
        ax_ret.plot(xs, ys, c=st.colors[5], zorder=3)
    plt.show()