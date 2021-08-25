import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad
import os
from random import shuffle, sample

from functions import *

if __name__ == "__main__":
    tau = 2.81
    j = 0.0574
    n_cha = 4
    n_clu = 10
    home = os.path.expanduser("~")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/"
    file = f"ca_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike =  f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    isis  = np.loadtxt(home + folder + file_spike)
    mean_isi = np.mean(isis)
    std_isi = np.std(isis)
    cv_isi = std_isi/mean_isi
    print(1/mean_isi, cv_isi)


    print(mean_isi, std_isi / mean_isi)
    spike_data_chunks = []
    chunks = []
    t = 0
    Tmax = 1_000
    for isi in isis:
        t += isi
        chunks.append(isi)
        if t > Tmax:
            spike_data_chunks.append(chunks.copy())
            chunks.clear()
            t = 0

    ws = np.logspace(-2, 1, 100)
    spectrum = []
    for w in ws:
        fws_real = []
        fws_img = []
        for isis in spike_data_chunks:
            fw = fourier_transformation_isis(w, isis)
            fws_real.append(fw.real)
            fws_img.append(fw.imag)
        spectrum.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))

    cas_theory = np.linspace(0, 1, 1001)
    dca = cas_theory[1] - cas_theory[0]
    p0s_theo_ca = []
    integral = 0
    for ca in reversed(cas_theory[1:]):
        print(ca)
        h = h_func(ca, tau, j, n_clu)
        d = d_func(ca, j, n_clu)
        if ca == 1:
            integral += 0
        elif ca >= 0.33:
            integral += np.exp(-h)*dca
        p0s_theo_ca.append(integral * np.exp(h) / d)

    norm = np.sum(p0s_theo_ca) * dca
    r0 = 1/norm

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    remove_top_right_axis([ax])
    ax.set_xlabel("$\omega$")
    ax.set_xscale("log")
    ax.set_ylabel("$S(\omega)$")
    ax.set_yscale("log")
    ax.axhline(r0, ls=":", c="C7", label="firing rate")
    ax.plot(ws, spectrum, c="C0", label="power spectrum")

    legend = ax.legend(loc=4, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    plt.savefig(home + "/Data/Calcium/Plots/10_ca_power_spectrum.pdf", transparent=True)
    plt.show()