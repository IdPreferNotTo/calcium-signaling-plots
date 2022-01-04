import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from random import shuffle, sample

from functions import *

if __name__ == "__main__":
    tau = 2.81
    j = 0.0728
    taua = 100
    ampa = 0.2
    n_cha = 4
    n_clu = 10
    home = os.path.expanduser("~")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/"
    file = f"ca_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    isis = np.loadtxt(home + folder + file_spike)
    mean_isi = np.mean(isis)
    std_isi = np.std(isis)
    cv_isi = std_isi / mean_isi
    print(1/mean_isi, cv_isi)

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    remove_top_right_axis([ax])

    k_correlations = []
    ks = np.arange(1, 6)
    var = k_corr(isis, isis, 0)
    for k in ks:
        covar  =k_corr(isis, isis, k)
        k_correlations.append(covar/var)

    shuffle(isis)
    k_correlations_shuffled = []
    var = k_corr(isis, isis, 0)
    for k in ks:
        covar_shuffled = k_corr(isis, isis, k)
        k_correlations_shuffled.append(covar_shuffled / var)

    ax.scatter(ks, k_correlations)
    ax.scatter(ks, k_correlations_shuffled)
    ax.axhline(0, ls=":", color="k")

    ax.set_xlabel("lag $k$")
    ax.set_xlim([0.8, 5.2])
    ax.set_ylabel(r"$\rho_k$")
    ax.set_ylim([-0.5, 0.2])

    plt.show()