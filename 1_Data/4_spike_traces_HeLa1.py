import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st

if __name__ == "__main__":
    nr = 1509
    histamine = 5.0
    home = os.path.expanduser("~")
    file_str = home + f"/Desktop/Ca data/Spikes/HeLa/1mircoMfuraTXT/{nr:d}_1_{histamine:.1f}histamine.dat"
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
        ax1.plot(times, data, lw=1, color=st.colors[4])
        plt.savefig(home + f"/Desktop/Ca data/Spikes/HeLa/1microMfuraPLT/HeLa_{nr:d}_1_{histamine:.1f}histamine_{i:d}.png")
        plt.show()
        plt.close()
