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

    for ax in axis:
        ax.set_ylim([0, 200])
        ax.set_xlabel("$i$")
        ax.set_ylabel("$T_i$")
    ax1.set_xlim([0, 25])
    ax2.set_xlim([0, 25])
    ax3.set_xlim([0, 50])
    ax4.set_xlim([0, 50])

    home = os.path.expanduser("~")
    with open(home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/spike_times.dat", "r") as infile:
        for idx, line in enumerate(infile):
            if idx == 27:
                ax = ax1
                y=20
            elif idx == 13:
                ax = ax2
                y = 20
            elif idx == 12:
                ax = ax3
                y = 40
            elif idx == 35:
                ax = ax4
                y = 40
            else:
                continue
            data_isi = line[:-2].split(" ")
            data_isi = [float(isi) for isi in data_isi]
            idxs = np.arange(len(data_isi))
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, data_isi, p0=(100, 150, 2))
            Tis = [popt[0] * np.exp(-idx / popt[2]) + popt[1] * (1 - np.exp(-idx / popt[2])) for idx in idxs]
            ax.plot(idxs, Tis, c="k", lw=1, label=rf"$n_{{\rm tr}}$ = {popt[2]:.1f}" + "\n" + rf"$\Delta T/T_0 = {(popt[1]-popt[0])/popt[0]:.1f}$")
            ax.axvspan(0, popt[2], facecolor="C7", alpha=0.5, zorder=0)
            ax.scatter(idxs, data_isi, fc="None", ec=st.colors[5], s=20)
            ax.axhline(popt[0], ls=":", c="C7")
            ax.axhline(popt[1], ls=":", c="C7")
            ax.annotate("", xy=(y, popt[1]), xytext=(y, popt[0]),
                        arrowprops=dict(arrowstyle='<->', connectionstyle="arc3", color="k", lw=1))
            ax.legend(fancybox=False, framealpha=1., edgecolor="k", fontsize=8)
    plt.show()