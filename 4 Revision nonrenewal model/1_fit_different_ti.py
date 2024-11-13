import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from scipy.stats import norm
from scipy.optimize import curve_fit

import functions as fc
import styles as st

def exponential_Ti_at_0(i, T0, T8, ntr):
    return T0*np.exp(-i/ntr) + T8*(1 - np.exp(-i/ntr))

def exponential_Ti_at_1(i, T1, T8, ntr):
    return T1*np.exp(-(i-1)/ntr) + T8*(1 - np.exp(-(i-1)/ntr))

def exponential_Ti_at_2(i, T2, T8, ntr):
    return T2*np.exp(-(i-2)/ntr) + T8*(1 - np.exp(-(i-2)/ntr))

if __name__ == "__main__":
    home = os.path.expanduser("~")
    idxs = [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 14, 16, 19, 20, 23, 27, 28, 29, 30, 21, 32, 33, 35]
    T0s = []
    dTs = []

    T0s_f1 = []
    T0s_f2 = []
    T0s_f3 = []
    err_T0s_f1 = []
    err_T0s_f2 = []
    err_T0s_f3 = []
    T8s_f1 = []
    T8s_f2 = []
    T8s_f3 = []
    dTs_f1 = []
    dTs_f2 = []
    dTs_f3 = []
    ntrs_f1 = []
    ntrs_f2 = []
    ntrs_f3 = []
    cvs_f1 = []
    cvs_f2 = []
    cvs_f3 = []

    f1 = open(home + "/Desktop/hek_fit/exp_statistics_t0.txt", "w")
    f2 = open(home + "/Desktop/hek_fit/exp_statistics_t0t1.txt", "w")
    f1.write("# idx, t0, cv, t8, ntr \n")
    f2.write("# idx, t0, cv, t8, ntr \n")
    for idx in idxs:
        Tis = np.loadtxt(home + f"/Data/calcium/experimental/Spikes/HEK/HEK2/spike_times_{idx:d}.dat")

        popt, pcov = curve_fit(exponential_Ti_at_0, [i for i in range(len(Tis))], Tis, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        perr = np.sqrt(np.diag(pcov))
        err_T0s_f1.append(perr[0])
        T0 = popt[0]
        T8 = popt[1]
        dT = T8 - T0
        ntr = popt[2]
        std_isi = np.std(Tis[int(2*popt[2]):])
        mean_isi = np.mean(Tis[int(2*popt[2]):])
        cv = std_isi/mean_isi
        T0s_f1.append(T0)
        T8s_f1.append(T8)
        dTs_f1.append(dT)
        ntrs_f1.append(ntr)
        cvs_f1.append(cv)

        popt, pcov = curve_fit(exponential_Ti_at_0, np.arange(1, len(Tis)), Tis[1:], bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        perr = np.sqrt(np.diag(pcov))
        err_T0s_f2.append(perr[0])
        T0 = popt[0]
        T8 = popt[1]
        dT = T8 - T0
        ntr = popt[2]
        std_isi = np.std(Tis[int(2*popt[2]):])
        mean_isi = np.mean(Tis[int(2*popt[2]):])
        cv = std_isi/mean_isi
        T0s_f2.append(T0)
        T8s_f2.append(T8)
        dTs_f2.append(dT)
        ntrs_f2.append(ntr)
        cvs_f2.append(cv)
        f1.write(f"{idx:d} {popt[0]:.1f} {std_isi / mean_isi:.2f} {popt[1]:.1f} {popt[2]:.2f} \n")

        st.set_default_plot_style()
        fig = plt.figure(tight_layout=True, figsize=(4, 2.5))
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])
        st.remove_top_right_axis([ax1])
        indexs = [i for i in range(len(Tis))]
        Tis_fit = [T0 * np.exp(-idx / ntr) + T8 * (1 - np.exp(-idx / ntr)) for idx in indexs]
        ax1.plot(indexs, Tis_fit, c="k", lw=1, label=f"$T_0$ = {popt[0]:.0f} $\pm$ {perr[0]:.0f}" + "\n" + f"$T_\infty$ = {popt[1]:.0f} $\pm$ {perr[1]:.0f}" + "\n" + rf"$n_{{\rm tr}}$ = {popt[2]:.1f} $\pm$ {perr[2]:.1f}" + "\n" + rf"$CV = {cv:.2f}$")
        ax1.set_xlabel("$i$")
        ax1.set_ylabel(r"$T_i$ / s")
        ax1.scatter(indexs, Tis, s=20, fc="w", ec=st.colors[7])
        ax1.axvspan(0, ntr, facecolor="C7", alpha=0.5, zorder=0)
        ax1.axhline(T0, ls=":", c="C7")
        ax1.axhline(T8, ls=":", c="C7")
        legend = ax1.legend(loc=4, fancybox=False, edgecolor="k", framealpha=0.0, prop={'size': 10})
        legend.get_frame().set_linewidth(0.5)
        plt.savefig(home + f"/Desktop/hek_fit/plots_no_t0/hek_fit_{idx:d}")
        plt.show()
        plt.close()

        popt, pcov = curve_fit(exponential_Ti_at_0, np.arange(2, len(Tis)), Tis[2:], bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        perr = np.sqrt(np.diag(pcov))
        err_T0s_f3.append(perr[0])
        T0 = popt[0] #popt[1] - (popt[1]-popt[0])*np.exp(2/popt[2])
        T8 = popt[1]
        dT = T8 - T0
        ntr = popt[2]
        std_isi = np.std(Tis[int(2*popt[2]):])
        mean_isi = np.mean(Tis[int(2*popt[2]):])
        cv = std_isi/mean_isi
        T0s_f3.append(T0)
        T8s_f3.append(T8)
        dTs_f3.append(dT)
        ntrs_f3.append(ntr)
        cvs_f3.append(cv)
        f2.write(f"{idx:d} {popt[0]:.1f} {std_isi / mean_isi:.2f} {popt[1]:.1f} {popt[2]:.2f} \n")

        st.set_default_plot_style()
        fig = plt.figure(tight_layout=True, figsize=(4, 2.5))
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])
        st.remove_top_right_axis([ax1])
        indexs = [i for i in range(len(Tis))]
        Tis_fit = [T0 * np.exp(-idx / ntr) + T8 * (1 - np.exp(-idx / ntr)) for idx in indexs]
        ax1.plot(indexs, Tis_fit, c="k", lw=1, label=f"$T_0$ = {popt[0]:.0f} $\pm$ {perr[0]:.0f}" + "\n" + f"$T_\infty$ = {popt[1]:.0f} $\pm$ {perr[1]:.0f}" + "\n" + rf"$n_{{\rm tr}}$ = {popt[2]:.1f} $\pm$ {perr[2]:.1f}" + "\n" + rf"$CV = {cv:.2f}$")
        ax1.set_xlabel("$i$")
        ax1.set_ylabel(r"$T_i$ / s")
        ax1.scatter(indexs, Tis, s=20, fc="w", ec=st.colors[7])
        ax1.axvspan(0, ntr, facecolor="C7", alpha=0.5, zorder=0)
        ax1.axhline(T0, ls=":", c="C7")
        ax1.axhline(T8, ls=":", c="C7")
        legend = ax1.legend(loc=4, fancybox=False, edgecolor="k", framealpha=0.0, prop={'size': 10})
        legend.get_frame().set_linewidth(0.5)
        plt.savefig(home + f"/Desktop/hek_fit/plots_no_t0t1/hek_fit_{idx:d}")
        plt.show()
        plt.close()

    width, height = st.figsize([1, 2], "large")
    fig, axes = plt.subplots(nrows=2, ncols=1, layout="constrained", figsize=(width, height))
    ax1 = axes[0]
    ax2 = axes[1]
    st.remove_top_right_axis([ax1, ax2])
    ax1.scatter(T0s_f1, err_T0s_f1, c=st.colors[1], ec="k", fc="w")
    for T0, errT0 in zip(T0s_f2, err_T0s_f2):
        xs = np.linspace(T0 - 3*errT0, T0 + 3*errT0, num= 1000)
        ys = norm.pdf(xs, T0, errT0)
    ax1.scatter(T0s_f2, err_T0s_f2, c=st.colors[1], ec="C2", fc="w")
    for T0, errT0 in zip(T0s_f3, err_T0s_f3):
        xs = np.linspace(T0 - 3*errT0, T0 + 3*errT0, num= 1000)
        ys = norm.pdf(xs, T0, errT0)
    ax2.scatter(T0s_f1, err_T0s_f1, c=st.colors[1], ec="k", fc="w")
    ax2.scatter(T0s_f2, err_T0s_f2, c=st.colors[1], ec="C3", fc="w")
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax1.set_xlabel("$T_0$")
    ax1.set_ylabel("Error $T_0$")
    plt.savefig(home + f"/Desktop/hek_fit/err_T0s.png", dpi=300)
    plt.show()
    plt.close()


    f1.close()
    f2.close()
    st.set_default_plot_style()
    width, height = st.figsize([2, 3], "large")
    fig, axes = plt.subplots(nrows=3, ncols=2, layout="constrained", figsize=(width, height))
    ax1 = axes[0, 0]
    ax2 = axes[0, 1]
    ax3 = axes[1, 0]
    ax4 = axes[1, 1]
    ax5 = axes[2, 0]
    ax6 = axes[2, 1]
    st.remove_top_right_axis([ax1, ax2, ax3, ax4, ax5, ax6])
    ax1.plot(T0s_f1, T0s_f1, c="k")
    ax1.scatter(T0s_f1, T0s_f2, s=20, ec="C2", fc="w", label="without $T_0$")
    ax1.set_ylabel("Estimate $T_0$")
    ax1.set_xlabel("$T_0$")
    ax1.set_ylim([0, 300])
    ax1.legend(fancybox=False, fontsize=8, loc=2)

    ax2.plot(T0s_f1, T0s_f1, c="k")
    ax2.scatter(T0s_f1, T0s_f3, s=20, ec="C3", fc="w", label="without $T_0, T_1$")
    ax2.set_ylabel("Estimate $T_0$")
    ax2.set_xlabel("$T_0$")
    ax2.set_ylim([0, 300])
    ax2.legend(fancybox=False, fontsize=8, loc=2)

    ax3.plot(dTs_f1, dTs_f1, c="k")
    ax3.scatter(dTs_f1, dTs_f2, s=20, ec="C2", fc="w", label="without $T_0$")
    ax3.set_ylabel("Estimate $\Delta T$")
    ax3.set_xlabel("$\Delta T$")
    ax3.set_ylim([0, 300])
    ax3.legend(fancybox=False, fontsize=8)

    ax4.plot(dTs_f1, dTs_f1, c="k")
    ax4.scatter(dTs_f1, dTs_f3, s=20, ec="C3", fc="w", label="without $T_0, T_1$")
    ax4.set_ylabel("Estimate $\Delta T$")
    ax4.set_xlabel("$\Delta T$")
    ax4.set_ylim([0, 300])
    ax4.legend(fancybox=False, fontsize=8)

    ax5.plot(ntrs_f1, ntrs_f1, c="k")
    ax5.scatter(ntrs_f1, ntrs_f2, s=20, ec="C2", fc="w", label="without $T_1$")
    ax5.set_ylabel(r"Estimate $n_{\rm tr}$")
    ax5.set_xlabel(r"$n_{\rm tr}$")
    ax5.set_ylim([0, 7])
    ax5.legend(fancybox=False, fontsize=8)


    ax6.plot(ntrs_f1, ntrs_f1, c="k")
    ax6.scatter(ntrs_f1, ntrs_f3, s=20, ec="C3", fc="w", label="without $T_0, T_1$")
    ax6.set_ylabel(r"Estimate $n_{\rm tr}$")
    ax6.set_xlabel(r"$n_{\rm tr}$")
    ax6.set_ylim([0, 7])
    ax6.legend(fancybox=False, fontsize=8)

    plt.savefig(home + f"/Desktop/hek_fit/estimate_parameters.png", dpi=300)
    plt.show()

