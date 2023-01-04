import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axis = [ax1, ax2, ax3, ax4]
    st.remove_top_right_axis(axis)

    home = os.path.expanduser("~")
    data = np.loadtxt(home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/HEK_transient_stationary_statistics.dat")
    idxs = data[:, 0]
    T0s = data[:, 1]
    T8s = data[:, 2]
    nTrs = data[:, 3]
    stdT8 = data[:, 4]
    varT0s = data[:, 5]
    varT8s = data[:, 6]
    varnTrs = data[:, 7]

    ax1.hist(T8s, color=st.colors[5], alpha=0.75)
    ax2.hist(stdT8/T8s, color=st.colors[5], alpha=0.75)
    ax3.hist(T8s - T0s, color=st.colors[5], alpha=0.75)
    ax4.hist(nTrs, color=st.colors[5], alpha=0.75)

    ax1.set_ylim([0, 11])
    ax2.set_ylim([0, 11])
    ax3.set_ylim([0, 11])
    ax4.set_ylim([0, 11])

    ax1.set_xlim([0, 400])
    ax2.set_xlim([0, 0.5])
    ax3.set_xlim([0, 400])
    ax4.set_xlim([0, 20])

    axin2 = ax2.inset_axes([0.5, 0.5, 0.45, 0.45])
    axin2.scatter(T8s, stdT8, s=20, facecolor="None", edgecolor=st.colors[5])
    popt, pcov = curve_fit(fc.linear_function, T8s, stdT8, p0=(-50, 0.2))
    print(popt)
    xs = np.linspace(0, 400)
    ys = popt[1] * (xs - popt[0])
    axin2.plot(xs, ys, c="k")
    axin2.set_xlim([0, 400])
    axin2.set_ylim([0, 100])
    axin2.set_xlabel(r"$T_\infty$")
    axin2.set_ylabel(r"$CV_T$")

    axin4 = ax4.inset_axes([0.5, 0.5, 0.45, 0.45])
    axin4.scatter(1/nTrs, T8s - T0s, s=20, facecolor="None", edgecolor=st.colors[5])
    print(nTrs)
    print((T8s - T0s)/T0s)
    #opt, pcov = curve_fit(fc.linear_function, T8s, stdT8, p0=(-50, 0.2))
    #print(popt)
    #xs = np.linspace(0, 400)
    #ys = popt[1] * (xs - popt[0])
    #axin2.plot(xs, ys, c="k")
    axin4.set_xlim([0, 2])
    axin4.set_ylim([0, 400])
    axin4.set_xlabel(r"$1/n_{\rm tr}$")
    axin4.set_ylabel(r"$\Delta T$")


    ax1.set_xlabel(r"$T_\infty$")
    ax2.set_xlabel(r"$CV_T$")
    ax3.set_xlabel(r"$\Delta T$")
    ax4.set_xlabel(r"$n_{\rm tr}$")

    ax1.set_ylabel(r"$p(T_\infty)$")
    ax2.set_ylabel(r"$p(CV_T)$")
    ax3.set_ylabel(r"$p(\Delta T)$")
    ax4.set_ylabel(r"$p(n_{\rm tr})$")

    data_small_small = []
    data_small_large = []
    data_large_small = []
    data_large_large = []
    for idx, T0, T8, nTr in zip(idxs, T0s, T8s, nTrs):
        dT  = (T8 - T0)/T0
        if nTr < 3 and dT < 2.:
            data_small_small.append(idx)
        elif nTr < 3 and dT > 2:
            data_small_large.append(idx)
        elif nTr > 3 and dT < 2:
            data_large_small.append(idx)
        else:
            data_large_large.append(idx)
    print(data_small_small)
    print(data_small_large)
    print(data_large_small)
    print(data_large_large)
    plt.show()
