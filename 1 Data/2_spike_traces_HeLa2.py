import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    nr = ["0609", "1409"][1]
    idxs = [1, 2]
    histamines = [[2, 5, 7.5, 10] ,[0.5, 1.0, 2.0, 5.0]][1]
    for idx in idxs:
        for histamine in histamines:
            home = os.path.expanduser("~")
            file_str = home + f"/Data/calcium/experimental/Spikes/HeLa/100microMHistamineTXT/{nr:s}_{idx:d}_{histamine:.1f}uM.dat"
            data = np.loadtxt(file_str)
            datas = np.transpose(data)
            times =  datas[0]

            for i, data in enumerate(datas[1:]):
                spiking: bool = False
                t_tmp: float = 0
                spike_times = []
                for t, ca in zip(times, data):
                    if t < 1000:
                        if ca > 0.15 and not spiking:
                            spike_times.append(t)
                            spiking = True
                        if ca < 0.15 and spiking:
                            spiking = False

                ISIs = []
                t0 = 0
                for t in spike_times:
                    ISIs.append(t - t0)
                    t0 = t

                st.set_default_plot_style()
                w = 3.25 * 1.25 / 2
                h = 2.66 * 1.25
                fig = plt.figure(tight_layout=True, figsize=(w, h))
                gs = gridspec.GridSpec(2, 1)
                ax1 = fig.add_subplot(gs[0])
                ax2 = fig.add_subplot(gs[1])
                st.remove_top_right_axis([ax1, ax2])
                ax1.text(0.6, 1.00, "HeLa cell", fontsize=11, transform=ax1.transAxes, ha="center", va='top')

                ax1.set_xlabel("$t$ / s")
                ax1.set_xlim([0, 1500])
                # ax1.axhline(0.4, lw=1, ls=":", color="C7")
                ax1.set_ylabel("$\Delta$F / a.u.")
                ax1.plot(times, data, lw=1, color=st.colors[7])

                nr_ISIs = len(ISIs)
                index_ISIs = np.arange(nr_ISIs)
                #popt, pcov = curve_fit(fc.exponential_Ti, index_ISIs, ISIs, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
                #i_ISI_fit = np.linspace(0, nr_ISIs)
                #ISI_fit = popt[0] * np.exp(-i_ISI_fit / popt[2]) + popt[1] * (1 - np.exp(-i_ISI_fit / popt[2]))

                #ax2.set_xlim([0, nr_ISIs])
                ax2.set_xlim([0, 10])
                ax2.set_xlabel("$i$")
                #ax2.set_ylim([0, 1.5 * popt[1]])
                ax2.set_ylabel("$T_i$ / s")
                ax2.scatter(index_ISIs[1:], ISIs[1:], fc="w", ec=st.colors[7], s=20, zorder=2)
                #ax2.plot(i_ISI_fit, ISI_fit, lw=1, c="k", zorder=3)
                #ax2.axvspan(0, popt[2], alpha=0.3, color="C7")
                #ax2.axhline(popt[0], ls=":", lw=1, c="k")
                #ax2.axhline(popt[1], ls=":", lw=1, c="k")


                plt.savefig(home + f"/Data/calcium/experimental/Spikes/HeLa/100microMHistaminePLT/HeLa_{nr:s}_{idx:d}_{histamine:.1f}histamine_{i:d}.png", dpi=300)
                plt.show()
                plt.close()
