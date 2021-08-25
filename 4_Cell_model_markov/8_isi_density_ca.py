import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
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
    home = os.path.expanduser("~")
    # Parameters
    tau = 2.81
    j = 0.0728 #0.0574
    n_cha = 4
    n_clu = 10
    print(tau, j)
    # Load Data
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/"
    file_spikes = "spike_times_markov_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu)
    isi_data = np.loadtxt(home + folder + file_spikes)

    mean_isi = np.mean(isi_data)
    std_isi = np.std(isi_data)
    cv_isi = std_isi/mean_isi
    cv2_isi = cv_isi**2
    ts = np.linspace(0, 50, 501)
    inv_gaus = []
    for t in ts:
        p = np.sqrt(mean_isi/(2*np.pi*cv2_isi*(t**3)))*np.exp(-(t - mean_isi)**2/(2*mean_isi*cv2_isi*t))
        inv_gaus.append(p)

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    #axins = inset_axes(ax, width="40%", height="40%", loc=1)
    remove_top_right_axis([ax])

    #axins.plot(cas_plot[:-1], fs_sim_ca, c="k", label="f(x)")
    #axins.plot(cas_plot, f_theo_ca, c="k", alpha=0.7, ls="--")
    #axins.plot(cas_plot[:-1], sqrt2ds_sim_ca, label="$\sqrt{2D(x)}$", c="C0", )
    #axins.plot(cas_plot, sqrt2D_theo_ca, c="C0", alpha=0.7, ls="--")
    #axins.legend()
    ax.plot(ts, inv_gaus, lw=1, c="k", label="Inverse Gaussian")
    ax.hist(isi_data, bins=50, density=True, alpha=.6, color="C0", label=f"$P(I)$ \n $\mu =$ {mean_isi:.1f} \n $C_V = {cv_isi:.2f}$")

    legend = ax.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    ax.set_ylabel(r"$P(I)$")
    ax.set_xlabel(r"$I$ [s]")
    plt.savefig(home + "/Data/Calcium/Plots/8_markov_isi_density.pdf", transparent=True)
    plt.show()
