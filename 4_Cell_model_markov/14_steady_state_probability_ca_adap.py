import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from functions import *

if __name__ == "__main__":
    # Parameters
    tau = 2.81
    j = 0.0728
    taua = 100
    ampa = 0.2
    home = os.path.expanduser("~")

    folder = "/Data/calcium_spikes_markov/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    isis = np.loadtxt(home + folder + file_spike)

    data = [set for set in data[1000:] if set[1] != 1]
    ts, cas, jpuffs, adaps = np.transpose(data)

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    remove_top_right_axis([ax])
    ax.set_ylabel(r"[Ca\textsuperscript{2+}]$_{\rm ER}$ [a.u.]")
    ax.set_xlabel(r"[Ca\textsuperscript{2+}]$_{\rm i}$ [a.u.]")
    ax.axvline(0.33, ls=":", lw=1, c="k")
    ax.axvline(1, ls="--", lw=1, c="k")
    ax.hist2d(cas, adaps, range = [[0.2, 1.05], [0.5, 1.05]], cmin = 0.01, bins=(50, 50), cmap="viridis", density=True)
    plt.show()