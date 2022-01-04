import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad
import os
from random import shuffle, sample

from functions import *

def mean_puff_single(x):
    m = 4
    n = 4
    r_opn = 0.13 * np.power(x / 0.33, 3) * ((1 + 0.33 ** 3) / (1 + x ** 3))
    r_ref = 1.3 * np.power(x / 0.33, 3) * ((1 + 0.33 ** 3) / (1 + x ** 3))
    r_cls = 50

    p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n)

    xs = [0, 0, 0, 0, 4, 3, 2, 1]
    mean = sum([x * p for x, p in zip(xs, p0s)])
    return mean


def means(N, tau, j):
    fs = []
    cas = np.linspace(0.01, 1.00, 100)
    for ca in cas:
        f = -(ca - 0.33)/tau + j*N*mean_puff_single(ca)
        fs.append(f)
    return fs


def intensity_puff_single(x):
    m = 4
    n = 4
    r_opn = 0.13 * np.power(x / 0.33, 3) * ((1. + 0.33 ** 3) / (1. + x ** 3))
    r_ref = 1.3 * np.power(x / 0.33, 3) * ((1. + 0.33 ** 3) / (1. + x ** 3))
    r_cls = 50

    xs = [0, 0, 0, 0, 4, 3, 2, 1]
    idxs = [0, 1, 2, 3, 4, 5, 6, 7]
    p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n)

    D_theory = 0
    for k in idxs:
        sum_over_i = 0
        f_from_k_to = f_from_k_invert_A(k, r_ref, r_opn, r_cls, m, n)
        for i in idxs:
            sum_over_i += xs[i] * f_from_k_to[i]
        D_theory += xs[k] * p0s[k] * sum_over_i
    return D_theory


def d_func(x, j, N):
    if x == 0:
        return 0
    else:
        return np.power(j, 2)*N*intensity_puff_single(x)


def f_func(x, tau, j, N):
    if x == 0:
        return -(x - 0.33) / tau
    else:
        return -(x - 0.33) / tau + j * N * mean_puff_single(x)


def g_func(x, tau, j, N):
    return f_func(x, tau, j, N)/d_func(x, j, N)


def h_func(x, tau, j, N):
    h = quad(g_func, 0.33, x, args=(tau, j, N))[0]
    return h


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
        for x in spike_data_chunks:
            fw = fourier_transformation_isis(w, x)
            fws_real.append(fw.real)
            fws_img.append(fw.imag)
        spectrum.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))

    shuffle(isis)
    spike_data_chunks_shuffle = []
    chunks_shuffle = []
    t = 0
    Tmax = 1_000
    for isi in isis:
        t += isi
        chunks_shuffle.append(isi)
        if t > Tmax:
            spike_data_chunks_shuffle.append(chunks_shuffle.copy())
            chunks_shuffle.clear()
            t = 0

    spectrum_shuffle = []
    for w in ws:
        fws_real = []
        fws_img = []
        for isis in spike_data_chunks_shuffle:
            fw = fourier_transformation_isis(w, isis)
            fws_real.append(fw.real)
            fws_img.append(fw.imag)
        spectrum_shuffle.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))


    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    remove_top_right_axis([ax])
    ax.set_xlabel("$\omega$")
    ax.set_xscale("log")
    ax.set_ylabel("$S(\omega)$")
    ax.set_yscale("log")
    ax.plot(ws, spectrum, c="C0", label="power spectrum")
    ax.plot(ws, spectrum_shuffle, c="C1", label="power spectrum")

    legend = ax.legend(loc=4, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    plt.savefig(home + "/Data/Calcium/Plots/12_ca_adap_power_spectrum.pdf", transparent=True)
    plt.show()