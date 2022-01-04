import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st

if __name__ == "__main__":
    nmax_channel = 5
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3, 4))
    gs = gridspec.GridSpec(nmax_channel, 1)
    axis = []
    for i in range(nmax_channel):
        ax = fig.add_subplot(gs[i])
        uni_data = np.random.randint(1, nmax_channel+1-i, 200)
        for k in range(1, nmax_channel+1-i):
            n = 0
            for x in uni_data:
                if x == k:
                    n+=1
            ax.bar(k, n/200, 0.8, align="center", color=st.colors[0])

        ax.set_xticklabels([])
        ax.set_xticks(range(1, nmax_channel+1))
        ax.set_xlim([0, nmax_channel+1])
        ax.text(1, 0.5, f"$n={nmax_channel-i:d}$", va="center", rotation="vertical", transform=ax.transAxes)
        axis.append(ax)

    st.remove_top_right_axis(axis)
    axis[int(nmax_channel/2)].set_ylabel(r"$P(n_{\rm opn})$")
    axis[-1].set_xticklabels(range(1, nmax_channel+1))
    axis[-1].set_xlabel(r"$n_{\rm opn}$")
    home = os.path.expanduser("~")
    plt.savefig(home + "/Data/Calcium/Plots/puff_gen_prob_nopen.pdf", transparent=True)
    plt.show()