import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st

if __name__ == "__main__":
    nr = [1509, 2309, 2709][0]
    idxs = [1, 2]
    histamines = [0.5, 1, 1.5, 2.0, 5.0]
    for idx in idxs:
        for histamine in histamines:
            home = os.path.expanduser("~")
            file_str = home + f"/Data/calcium/experimental/Spikes/HeLa/1mircoMfuraTXT/{nr:d}_{idx:d}_{histamine:.1f}histamine.dat"
            data = np.loadtxt(file_str)
            datas = np.transpose(data)
            times =  datas[0]

            for i, data in enumerate(datas[1:]):
                st.set_default_plot_style()
                fig = plt.figure(tight_layout=True, figsize=(4, 2))
                gs = gridspec.GridSpec(1, 1)
                ax1 = fig.add_subplot(gs[0])
                st.remove_top_right_axis([ax1])

                ax1.set_xlabel("$t$ / s")
                ax1.set_ylabel("Ratio (340/380)")
                ax1.set_xlim([0, 1500])
                ax1.plot(times, data, lw=1, color=st.colors[7])
                plt.savefig(home + f"/Data/calcium/experimental/Spikes/HeLa/1microMfuraPLT/HeLa_{nr:d}_{idx:d}_{histamine:.1f}histamine_{i:d}.png", dpi=300, transparent=True)
                plt.show()
                plt.close()
