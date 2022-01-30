import numpy as np
from random import shuffle, sample
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from functions import *


if __name__ == "__main__":
    home = os.path.expanduser("~")
    tau = 2.81
    j = 0.489
    taua = 100
    ampa = 0.2

    spike_data = []
    for run in range(10):
        home = os.path.expanduser("~")

        folder = "/CLionProjects/PhD/calcium_spikes_langevin/out/"
        file = f"ca_langevin_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
        file_spike = f"spike_times_langevin_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
        data = np.loadtxt(home + folder + file_spike)
        spike_data.append(data)

    spike_data = [isi for data in spike_data for isi in data]
    # Cut transient
    spike_data = spike_data[100:]

    sum_corr = 0
    var = k_corr(spike_data, spike_data, 0)
    k_covar = []
    ks = np.arange(1, 5)
    for k in ks:
        covar = k_corr(spike_data, spike_data, k)
        k_covar.append(covar/var)
        sum_corr += covar/var

    mean_isi = np.mean(spike_data)
    std_isi = np.std(spike_data)

    print(mean_isi, std_isi/mean_isi)
    spike_data_chunks = []
    chunks = []
    t = 0
    Tmax = 10_000
    for isi in spike_data:
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
        spectrum.append((1./Tmax) * (np.var(fws_real) + np.var(fws_img)))

    shuffle(spike_data)

    sum_cor_shuffledr = 0
    var = k_corr(spike_data, spike_data, 0)
    k_covar_shuffled = []
    ks = np.arange(1, 5)
    for k in ks:
        covar_shuffled = k_corr(spike_data, spike_data, k)
        k_covar_shuffled.append(covar/var)
        sum_corr += covar/var

    spike_data_shuffled_chunks = []
    chunks_shuffled = []
    t = 0
    Tmax = 10_000
    for isi in spike_data:
        t += isi
        chunks_shuffled.append(isi)
        if t > Tmax:
            spike_data_shuffled_chunks.append(chunks_shuffled.copy())
            chunks_shuffled.clear()
            t = 0

    spectrum_shuffled = []
    for w in ws:
        fws_real = []
        fws_img = []
        for isis in spike_data_shuffled_chunks:
            fw = fourier_transformation_isis(w, isis)
            fws_real.append(fw.real)
            fws_img.append(fw.imag)
        spectrum_shuffled.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))

    fig, ax = plt.subplots(1,1)
    ax.plot(ws, spectrum, c="k", label="unshuffled")
    ax.plot(ws, spectrum_shuffled, c="C7", label="shuffled")
    ax.set_xlabel("$\omega$")
    ax.set_ylabel("$S(\omega)$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(prop={"size": 7}, loc=1, ncol=1, framealpha=1., edgecolor="k")
    leg = ax.get_legend()
    leg.get_frame().set_linewidth(0.5)
    leg.legendHandles[0].set_color('k')
    leg.legendHandles[1].set_color('k')


    axins = inset_axes(ax, width="40%", height="30%", loc='lower left',
                   bbox_to_anchor=(0.55,0.1,1,1), bbox_transform=ax.transAxes)
    axins.scatter(ks, k_covar, c="k", zorder=4)
    axins.scatter(ks, k_covar_shuffled, c="C7", zorder=3)
    axins.set_xlabel("$k$")
    axins.set_ylabel(r"$\rho_k$")
    axins.axhline(0, ls=":", c="k")
    axins.set_ylim([-0.3, 0.1])
    plt.show()



