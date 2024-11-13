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
    # fig = plt.figure(tight_layout=True, figsize=(4, 3))
    # gs = gridspec.GridSpec(2, 2)
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax2 = fig.add_subplot(gs[0, 1])
    # ax3 = fig.add_subplot(gs[1, 0])
    # ax4 = fig.add_subplot(gs[1, 1])
    w = 3.25 * 1.25
    h = 3.75 * 1.25
    fig, axs = plt.subplots(nrows=3, ncols=2, layout="constrained", figsize=(w, h))
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[1, 0]
    ax4 = axs[1, 1]
    ax5 = axs[2, 0]
    ax6 = axs[2, 1]

    st.remove_top_right_axis([ax1, ax2, ax3, ax4, ax5, ax6])
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"C", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"D", fontsize=11, transform=ax4.transAxes, va='top')
    ax5.text(0.05, 0.95, r"E", fontsize=11, transform=ax5.transAxes, va='top')
    ax6.text(0.05, 0.95, r"F", fontsize=11, transform=ax6.transAxes, va='top')

    ax1.set_ylabel(r"$T_i$ / s")
    ax1.set_xlabel(r"$i$")
    ax2.set_ylabel(r"$T_i$ / s")
    ax2.set_xlabel(r"$i$")
    ax3.set_ylabel(r"$p(n_{\rm tr})$")
    ax3.set_xlabel(r"$n_{\rm tr}$")
    ax4.set_ylabel(r"$p(\Delta T)$")
    ax4.set_xlabel(r"$\Delta T$  / s")
    ax5.set_ylabel(r"$p(\tau_{\rm er})$")
    ax5.set_xlabel(r"$\tau_{\rm er}$ / min")
    ax6.set_ylabel(r"$p(\varepsilon)$")
    ax6.set_xlabel(r"$\varepsilon$")

    #------------------ EXAMPLES ---------------------------------------------------------------------------------------
    home = os.path.expanduser("~")
    ax1.set_xlim([0.5, 15])
    ax1.set_xticks([i for i in range(1, 16, 2)])
    ax1.set_ylim([0, 400])
    ax2.set_xlim([0.5, 30])
    ax2.set_xticks([i for i in range(1, 31, 4)])
    ax2.set_ylim([0, 400])
    for ax, cell_idx in zip([ax1, ax2], [9, 11]):
        exp_isi = np.loadtxt(home + f"/Data/calcium/experimental/Spikes/HEK/HEK2/spike_times_{cell_idx:d}.dat")
        exp_idx = [i for i in range(len(exp_isi))]
        popt, pcov = curve_fit(fc.exponential_Ti, exp_idx, exp_isi, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        T0 = popt[0]
        T8 = popt[1]
        n_tr = popt[2]
        mean_exp_isi = np.mean(exp_isi[int(2 * popt[2]):])
        std_exp_isi = np.std(exp_isi[int(2 * popt[2]):])
        cv_exp_isi = std_exp_isi / mean_exp_isi
        print(cv_exp_isi)
        Tis = [popt[0] * np.exp(-idx / popt[2]) + popt[1] * (1 - np.exp(-idx / popt[2])) for idx in exp_idx]

        ax.plot(exp_idx, Tis, c="k", lw=1, label=rf"$n_{{\rm tr}} = {popt[2]:.1f}$" + "\n"
                                                 + rf"$\Delta T = {round(popt[1] - popt[0],-1):.0f}$s")
        ax.scatter(exp_idx[1:], exp_isi[1:], s=15, fc="w", ec=st.colors[7])
        ax.axhline(popt[0], ls=":", c="C7")
        ax.axhline(popt[1], ls=":", c="C7")
        ax.axvspan(0, math.ceil(2*popt[2]), facecolor="C7", alpha=0.3, zorder=0)
        leg = ax.legend(fancybox=False, framealpha=1., edgecolor="k", fontsize=7)
        leg.get_frame().set_linewidth(0.75)
        for t in leg.get_texts():
            t.set_ha("center")

    # ----------------- STATISTICS -------------------------------------------------------------------------------------
    statistics = np.asarray([[40, 0.12, 1.1, 343],
         [24, 0.11, 4.6, 121],
         [21, 0.12, 7.4, 224],
         [13, 0.18, 2.1, 219],
         [25, 0.20, 3.5, 187],
         [19, 0.12, 2.4, 284],
         [29, 0.31, 0.9, 187],
         [29, 0.17, 2.4, 311],
         [29, 0.17, 1.6, 297],
         [29, 0.17, 6.4, 146],
         [26, 0.07, 3.3, 139],
         [21, 0.15, 3.8, 136],
         [22, 0.12, 4.9, 164],
         [16, 0.27, 1.8, 123],
         [35, 0.11, 1.6, 91],
         [20, 0.26, 1.8, 165],
         [21, 0.20, 3.8, 403],
         [41, 0.1, 2.6, 143],
         [46, 0.09, 4.3, 194],
         [46, 0.09, 4.2, 194],
         [31, 0.14, 1.9, 174],
         [36, 0.10, 7.4, 232],
         [36, 0.14, 4.7, 109],
         [15, 0.15, 5.8, 175]])

    parameters_it0 = np.asarray([[10.2, 8.82e-03, 4.14e+02, 1.78e-01],
                 [14.1, 9.90e-03, 6.49e+02, 5.78e-02],
                 [11.6, 1.16e-02, 1.94e+03, 4.53e-02],
                 [5.34, 2.14e-02, 5.12e+02, 1.44e-01],
                 [4.20, 1.86e-02, 6.50e+02, 3.68e-02],
                 [12.9, 1.19e-02, 8.02e+02, 1.60e-01],
                 [2.71, 2.53e-02, 1.75e+02, 8.52e-02],
                 [5.48, 1.47e-02, 7.59e+02, 6.30e-02],
                 [5.48, 1.47e-02, 5.01e+02, 9.45e-02],
                 [5.48, 1.47e-02, 9.81e+02, 1.90e-02],
                 [52.1, 6.80e-03, 6.00e+02, 1.15e-01],
                 [6.46, 1.51e-02, 5.79e+02, 5.23e-02],
                 [11.4, 1.13e-02, 9.36e+02, 5.84e-02],
                 [2.43, 3.13e-02, 2.48e+02, 6.63e-02],
                 [12.2, 8.48e-03, 1.54e+02, 8.71e-02],
                 [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
                 [3.93, 2.05e-02, 1.60e+03, 4.10e-02],
                 [10.3, 8.65e-03, 4.15e+02, 5.27e-02],
                 [18.6, 6.02e-03, 9.31e+02, 4.84e-02],
                 [18.6, 6.02e-03, 9.40e+02, 4.83e-02],
                 [7.56, 1.17e-02, 3.65e+02, 8.35e-02],
                 [15.5, 7.47e-03, 1.92e+03, 3.41e-02],
                 [7.71, 1.09e-02, 5.37e+02, 2.38e-02],
                 [6.94, 1.74e-02, 1.19e+03, 5.05e-02]])
    parameters_it2 = np.asarray([[6.09, 1.24e-02, 4.20e+02, 1.09e-01],
         [19.5, 8.95e-03, 7.06e+02, 6.44e-02],
         [22.4, 9.54e-03, 2.10e+03, 6.01e-02],
         [3.22, 2.73e-02, 5.18e+02, 9.78e-02],
         [7.11, 1.32e-02, 7.18e+02, 5.56e-02],
         [9.65, 1.31e-02, 8.00e+02, 1.31e-01],
         [1.84, 3.53e-02, 1.83e+02, 6.42e-02],
         [5.95, 1.40e-02, 7.87e+02, 6.71e-02],
         [4.42, 1.72e-02, 4.83e+02, 8.02e-02],
         [9.92, 1.04e-02, 1.04e+03, 3.13e-02],
         [53.3, 6.89e-03, 6.04e+02, 1.17e-01],
         [9.39, 1.26e-02, 5.91e+02, 6.72e-02],
         [23.5, 9.06e-03, 1.03e+03, 8.33e-02],
         [2.41, 3.14e-02, 2.28e+02, 6.75e-02],
         [10.2, 9.23e-03, 1.56e+02, 7.55e-02],
         [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
         [7.24, 1.42e-02, 1.69e+03, 6.79e-02],
         [13.9, 7.35e-03, 3.82e+02, 6.55e-02],
         [28.5, 5.03e-03, 9.54e+02, 5.86e-02],
         [35.4, 4.72e-03, 9.65e+02, 6.54e-02],
         [8.56, 1.08e-02, 3.53e+02, 9.42e-02],
         [32.3, 5.80e-03, 2.09e+03, 4.88e-02],
         [12.4, 8.30e-03, 5.46e+02, 3.36e-02],
         [77.0, 1.11e-02, 1.55e+03, 1.09e-01]])

    ntrs = [s[2] for s in statistics]
    dTs = [s[3] for s in statistics]
    taus = np.asarray([p[0] for p in parameters_it2])
    ps = np.asarray([p[1] for p in parameters_it2])
    tau_ers = np.asarray([p[2] for p in parameters_it2])
    eps_ers = np.asarray([p[3] for p in parameters_it2])



    # ax1.set_xlim([100, 3000])
    # ax1.set_xscale("log")
    ax3.set_yticklabels([])
    ax3.hist(ntrs, density=True, bins = [0, 2, 4, 6, 8, 10], color=st.colors[8])
    ax3.axvline(np.mean(ntrs), c="k", label=r"$\mu(\tau_{\rm er})$")
    ax3.axvline(np.mean(ntrs) - np.std(ntrs), c="k", ls=":", label=r"$\mu(\tau_{\rm er}) \pm \sigma(\tau_{\rm er})$")
    ax3.axvline(np.mean(ntrs) + np.std(ntrs), c="k", ls=":")

    ax4.set_yticklabels([])
    ax4.hist(dTs, density=True, bins=[0, 100, 200, 300, 400, 500], color=st.colors[8])
    ax4.axvline(np.mean(dTs), c="k", label=r"$\mu(\tau_{\rm er})$")
    ax4.axvline(np.mean(dTs) - np.std(dTs), c="k", ls=":", label=r"$\mu(\tau_{\rm er}) \pm \sigma(\tau_{\rm er})$")
    ax4.axvline(np.mean(dTs) + np.std(dTs), c="k", ls=":")

    ax5.set_xlim([0, 30])
    ax5.set_yticklabels([])
    ax5.hist(tau_ers/60, density=True, bins=[0, 5, 10, 15, 20, 25, 30], color=st.colors[8])
    ax5.axvline(np.mean(tau_ers/60), c="k", label=r"$\mu(\tau_{\rm er})$")
    ax5.axvline(np.mean(tau_ers/60) - np.std(tau_ers/60), c="k", ls=":", label=r"$\mu(\tau_{\rm er}) \pm \sigma(\tau_{\rm er})$")
    ax5.axvline(np.mean(tau_ers/60) + np.std(tau_ers/60), c="k", ls=":")

    ax6.set_xlim([0, 0.2])
    ax6.set_yticklabels([])
    ax6.hist(eps_ers, density=True, bins=[0, 0.04, 0.08, 0.12, 0.16, 0.2], color=st.colors[8])
    ax6.axvline(np.mean(eps_ers), c="k", label=r"$\mu(\varepsilon)$")
    ax6.axvline(np.mean(eps_ers) - np.std(eps_ers), c="k", ls=":", label=r"$\mu(\varepsilon) \pm \sigma(\varepsilon)$")
    ax6.axvline(np.mean(eps_ers) + np.std(eps_ers), c="k", ls=":")
    # ax4.legend(fancybox=False, fontsize=8)

    print("n_tr: ",np.mean(ntrs), np.std(ntrs, ddof=1), np.std(ntrs, ddof=1)/np.mean(ntrs))
    print("dT: ",np.mean(dTs), np.std(dTs, ddof=1), np.std(dTs, ddof=1)/np.mean(dTs))
    print("tau_er: ",np.mean(tau_ers), np.std(tau_ers, ddof=1), np.std(tau_ers, ddof=1)/np.mean(tau_ers))
    print("e_er: ",np.mean(eps_ers), np.std(eps_ers, ddof=1), np.std(eps_ers, ddof=1)/np.mean(eps_ers))
    #plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/SUB2/figures/fig11.pdf", dpi=300, transparent=True)

    plt.show()