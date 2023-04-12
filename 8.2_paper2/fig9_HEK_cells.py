import os

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
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"C", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"D", fontsize=11, transform=ax4.transAxes, va='top')

    ax1.set_ylabel(r"$T_i$")
    ax1.set_xlabel(r"$i$")
    ax2.set_ylabel(r"$T_i$")
    ax2.set_xlabel(r"$i$")
    ax3.set_ylabel(r"$p(\tau_{\rm er})$")
    ax3.set_xlabel(r"$\tau_{\rm er}$")
    ax4.set_ylabel(r"$p(\varepsilon)$")
    ax4.set_xlabel(r"$\varepsilon$")

    #------------------ EXAMPLES ---------------------------------------------------------------------------------------
    home = os.path.expanduser("~")
    ax1.set_xlim([0, 15])
    ax1.set_ylim([0, 400])
    ax2.set_xlim([0, 30])
    ax2.set_ylim([0, 400])
    for ax, cell_idx in zip([ax1, ax2], [9, 11]):
        exp_isi = np.loadtxt(home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/spike_times_{cell_idx:d}.dat")
        exp_idx = [i for i in range(len(exp_isi))]
        t0 = exp_isi[0]
        func = lambda x, t8, ntr: fc.exponential_Ti(x, t0, t8, ntr)
        popt, pcov = curve_fit(func, exp_idx, exp_isi, p0=(exp_isi[-1], 5.))
        perr = np.sqrt(np.diag(pcov))
        popt = np.insert(popt, 0, t0)
        perr = np.insert(perr, 0, 0)
        mean_exp_isi = np.mean(exp_isi[int(2 * popt[2]):])
        std_exp_isi = np.std(exp_isi[int(2 * popt[2]):])
        cv_exp_isi = std_exp_isi / mean_exp_isi
        Tis = [popt[0] * np.exp(-idx / popt[2]) + popt[1] * (1 - np.exp(-idx / popt[2])) for idx in exp_idx]

        ax.plot(exp_idx, Tis, c="k", lw=1, label=rf"$n_{{\rm tr}}={popt[2]:.1f}$" + "\n" + rf"$\Delta T = {popt[1]-popt[0]:.0f}$")
        ax.scatter(exp_idx, exp_isi, s=15, fc="w", ec=st.colors[5])
        ax.axhline(popt[0], ls=":", c="C7")
        ax.axhline(popt[1], ls=":", c="C7")
        ax.axvspan(0, popt[2], facecolor="C7", alpha=0.5, zorder=0)
        leg = ax.legend(fancybox=False, framealpha=1., edgecolor="k", fontsize=7)
        leg.get_frame().set_linewidth(0.75)
        for t in leg.get_texts():
            t.set_ha("right")

    # ----------------- STATISTICS -------------------------------------------------------------------------------------
    parameter2 = [[1.01e+01, 8.46e-03, 4.79e+02, 2.83e-01],
                  [1.25e+01, 9.07e-03, 6.69e+02, 1.05e-01],
                  [7.50e+00, 1.19e-02, 1.41e+03, 5.84e-02],
                  [5.42e+00, 1.41e-02, 5.06e+02, 9.97e-02],
                  [4.82e+00, 1.54e-02, 6.88e+02, 5.27e-02],
                  [1.02e+01, 9.32e-03, 8.46e+02, 1.68e-01],
                  [3.04e+00, 2.24e-02, 1.91e+02, 1.31e-01],
                  [5.74e+00, 1.35e-02, 7.59e+02, 1.06e-01],
                  [5.74e+00, 1.35e-02, 5.03e+02, 1.60e-01],
                  [6.68e+00, 1.17e-02, 8.27e+02, 3.53e-02],
                  [7.32e+01, 5.05e-03, 9.10e+02, 2.43e-01],
                  [6.90e+00, 1.19e-02, 5.70e+02, 6.29e-02],
                  [8.72e+00, 1.02e-02, 8.46e+02, 6.76e-02],
                  [3.45e+00, 2.01e-02, 2.38e+02, 6.10e-02],
                  [1.20e+01, 8.53e-03, 1.76e+02, 2.01e-01],
                  [3.60e+00, 1.94e-02, 3.04e+02, 7.38e-02],
                  [5.66e+00, 1.32e-02, 1.61e+03, 5.26e-02],
                  [1.01e+01, 8.46e-03, 4.04e+02, 1.00e-01],
                  [1.97e+01, 6.36e-03, 9.73e+02, 1.48e-01],
                  [1.97e+01, 6.36e-03, 9.46e+02, 1.50e-01],
                  [7.94e+00, 1.03e-02, 3.53e+02, 1.41e-01],
                  [1.57e+01, 7.55e-03, 1.63e+03, 9.84e-02],
                  [7.94e+00, 1.03e-02, 5.16e+02, 4.67e-02],
                  [7.19e+00, 1.11e-02, 1.03e+03, 3.96e-02]]
    tau_ers = [p[2] for p in parameter2]
    eps_ers = [p[3] for p in parameter2]

    ax3.set_xlim([0, 2000])
    # ax1.set_xlim([100, 3000])
    # ax1.set_xscale("log")
    ax3.set_yticklabels([])
    ax3.hist(tau_ers, density=True, bins=3, color=st.colors[4])
    ax3.axvline(np.mean(tau_ers), c="k", label=r"$\mu(\tau_{\rm er})$")
    ax3.axvline(np.mean(tau_ers) - np.std(tau_ers), c="k", ls=":", label=r"$\mu(\tau_{\rm er}) \pm \sigma(\tau_{\rm er})$")
    ax3.axvline(np.mean(tau_ers) + np.std(tau_ers), c="k", ls=":")
    # ax1.legend(fancybox=False, fontsize=8)


    ax4.set_xlim([0, 1])
    # ax2.set_xlim([0.01, 1])
    # ax2.set_xscale("log")
    ax4.set_yticklabels([])
    ax4.hist(eps_ers, density=True, bins=3, color=st.colors[4])
    ax4.axvline(np.mean(eps_ers), c="k", label=r"$\mu(\varepsilon)$")
    ax4.axvline(np.mean(eps_ers) - np.std(eps_ers), c="k", ls=":", label=r"$\mu(\varepsilon) \pm \sigma(\varepsilon)$")
    ax4.axvline(np.mean(eps_ers) + np.std(eps_ers), c="k", ls=":")
    # ax4.legend(fancybox=False, fontsize=8)

    print(np.mean(tau_ers), np.std(tau_ers)/np.mean(tau_ers), np.mean(eps_ers), np.std(eps_ers)/np.mean(eps_ers))
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig9.pdf", transparent=True)

    plt.show()