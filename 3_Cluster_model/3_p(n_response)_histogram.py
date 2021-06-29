import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

if __name__ == "__main__":
    fig = plt.figure(tight_layout=True, figsize=(9/2, 6))
    gs = gridspec.GridSpec(10, 1)
    axis = []
    for i in range(10):
        ax = fig.add_subplot(gs[i])
        uni_data = np.random.randint(1, 11-i, 100)
        for k in range(1, 11-i):
            n = 0
            for x in uni_data:
                if x == k:
                    n+=1
            ax.bar(k, n/100, 0.8, align="center", color="C0")

        ax.set_xticklabels([])
        ax.set_xticks([1,2,3,4,5,6,7,8,9,10])
        ax.set_xlim([0, 11])
        ax.yaxis.set_label_position("right")
        ax.set_ylabel("n={:d}".format(10-i))
        axis.append(ax)
    axis[-1].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
    home = os.path.expanduser("~")
    plt.savefig(
        home + "/Data/Calcium/Plots/histo_responding_channel.pdf", transparent=True)
    plt.show()