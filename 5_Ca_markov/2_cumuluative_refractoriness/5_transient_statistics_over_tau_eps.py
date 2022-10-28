import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axis_top = [ax1, ax2]
    axis_bot = [ax3, ax4]
    axis = axis_top + axis_bot
    st.remove_top_right_axis(axis)

    ax1.set_ylabel(r"$n_{\rm tr}$")
    ax1.set_xlabel(r"$\tau_{er}$")
    ax1.set_xscale("log")

    ax2.set_ylabel(r"$\Delta T$")
    ax2.set_xlabel(r"$\tau_{er}$")
    ax2.set_xscale("log")

    ax3.set_ylabel(r"$n_{\rm tr}$")
    ax3.set_xlabel(r"$\varepsilon_{er}$")
    ax3.set_xscale("log")

    ax4.set_ylabel(r"$\Delta T$")
    ax4.set_xlabel(r"$\varepsilon_{er}$")
    ax4.set_xscale("log")
    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.1
    tau_er_fix = 100
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)

    for tau, p, in zip(taus, ps):
        dTs = []
        n_trs = []
        n_trs_theory = []
        for tau_er in tau_ers:
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er_fix)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=[100, 150, 2])
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs.append(T8 - T0)
            n_trs.append(n_tr)
            n_trs_theory.append(-1/np.log((1-eps_er_fix)*np.exp(-T8/tau_er)))
        ax1.scatter(tau_ers, n_trs, fc="w", ec=st.colors[0], alpha=0.5, s=20, zorder=3)
        ax2.scatter(tau_ers, dTs, fc="w", ec=st.colors[0], alpha=0.5, s=20, zorder=3)
        ax1.plot(tau_ers, n_trs_theory, c=st.colors[5])

        dTs = []
        n_trs = []
        n_trs_theory = []
        dTs_theory = []
        for eps_er in eps_ers:
            print(eps_er)
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er_fix, eps_er)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=(100, 150, 2))
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs.append(T8 - T0)
            n_trs.append(n_tr)
            n_trs_theory.append(-1 / np.log((1 - eps_er) * np.exp(-T8 / tau_er_fix)))
        ax3.scatter(eps_ers, n_trs, fc="w", ec=st.colors[0], alpha=0.5, s=20, zorder=3)
        ax4.scatter(eps_ers, dTs, fc="w", ec=st.colors[0], alpha=0.5, s=20, zorder=3)

        ax3.plot(eps_ers, n_trs_theory, c=st.colors[5])
    plt.show()