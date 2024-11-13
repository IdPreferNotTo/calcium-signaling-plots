import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors

import styles as st
import functions as fc

def tau_of_p(x, a, b):
    return b/(x-a)

def p_of_tau(x, a, b):
    return a + b/x

if __name__ == "__main__":
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 1.55 * 1.25
    fig, axes = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(w, h))
    ax1 = axes[0]
    ax2 = axes[1]
    st.remove_top_right_axis([ax1, ax2])
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    ax1.set_ylim([0, 0.05])
    ax2.set_ylim([0, 0.05])

    statistics = np.asarray([
        [40, 0.12, 1.1, 343],
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

    parameters_it2 = np.asarray(
        [[6.09, 1.24e-02, 4.20e+02, 1.09e-01],
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

    pis = parameters_it2[:, 1]
    tauis = parameters_it2[:, 0]
    tauers = parameters_it2[:, 2]
    epsers = parameters_it2[:, 3]

    t0s = statistics[:, 0]
    cvs = statistics[:, 1]
    t8s = statistics[:, 3]
    ntrs = statistics[:, 2]


    cerfix_0 = 1 / (1 + epsers / ( 1- epsers/2) * tauers / t8s)
    # print(cerfix_0)
    print(cerfix_0)

    #ax1.set_xscale("log")
    #ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\tau$ / s")
    ax1.set_ylabel(r"$p$")
    ax1.set_xlim([1, 100])
    ax1.set_ylim([0.001, 0.1])
    ax1.scatter(tauis, pis, fc="w", s=20, ec=st.colors[8], zorder=3)

    #ax1.set_xscale("log")
    #ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$\tau$ / s")
    ax2.set_ylabel(r"$p \langle c_{\rm er}^* \rangle$")
    ax2.set_xlim([1, 100])
    ax2.set_ylim([0.001, 0.1])
    ax2.scatter(tauis, cerfix_0*pis, fc="w", s=20, ec=st.colors[8], zorder=3)


    popt, pcov = curve_fit(p_of_tau, tauis, pis)
    a = popt[0]
    b = popt[1]
    da = pcov[0,0]
    db = pcov[1,1]
    xs = np.linspace(2, 50, num=19)
    ys = p_of_tau(xs, a, b)
    ax1.plot(xs, ys, c="k", ls=":", zorder=2, label=r"$p(\tau) = a + b/\tau$")

    popt, pcov = curve_fit(p_of_tau, tauis, cerfix_0*pis)
    a = popt[0]
    b = popt[1]
    da = pcov[0,0]
    db = pcov[1,1]
    print(a, b, da, db)
    ys = p_of_tau(xs, a, b)
    ax2.plot(xs, ys, c="k", ls=":", zorder=2, label=r"$p(\tau) = a + b/\tau$")

    mu = fc.mean_jp_single_theory(0.5, 5, 3, 1)
    ys2  = (0.5-0.2)/(10*mu*xs)
    ax1.plot(xs, ys2, c="k", zorder=2, label=r"$p_{\rm bif}(\tau)$")
    ax2.plot(xs, ys2, c="k", zorder=2, label=r"$p_{\rm bif}(\tau)$")

    taus_bif = xs
    ps_bif = ys2
    print(taus_bif)
    print(ps_bif)


    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.55, 0.0), loc=3,
                     ncol=2, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(0.75)

    plt.savefig(home +"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/SUB2/figures/fig15.pdf")
    plt.show()
