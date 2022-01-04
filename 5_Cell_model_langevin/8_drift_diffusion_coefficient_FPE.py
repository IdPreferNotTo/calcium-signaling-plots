import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import os
from functions import *


def steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n):
    A = np.array([[-n * r_ref, 0, 0, 0, 0, 0, 0, r_cls],
                  [n * r_ref, -n * r_ref, 0, 0, 0, 0, 0, 0],
                  [0, n * r_ref, -n * r_ref, 0, 0, 0, 0, 0],
                  [0, 0, n * r_ref, -n * r_opn, 0, 0, 0, 0],
                  [0, 0, 0, r_opn, -r_cls, 0, 0, 0],
                  [0, 0, 0, r_opn, r_cls, -r_cls, 0, 0],
                  [0, 0, 0, r_opn, 0, r_cls, -r_cls, 0],
                  [1, 1, 1, 1, 1, 1, 1, 1]])
    Ainv = np.linalg.inv(A)
    inhomgeneity = np.array([0, 0, 0, 0, 0, 0, 0, 1])
    p0s = Ainv.dot(inhomgeneity)
    return p0s


def f_from_k_invert_A(k, r_ref, r_opn, r_cls, m, n):
    A = np.array([[-n * r_ref, 0, 0, 0, 0, 0, 0, r_cls],
                  [n * r_ref, -n * r_ref, 0, 0, 0, 0, 0, 0],
                  [0, n * r_ref, -n * r_ref, 0, 0, 0, 0, 0],
                  [0, 0, n * r_ref, -n * r_opn, 0, 0, 0, 0],
                  [0, 0, 0, r_opn, -r_cls, 0, 0, 0],
                  [0, 0, 0, r_opn, r_cls, -r_cls, 0, 0],
                  [0, 0, 0, r_opn, 0, r_cls, -r_cls, 0],
                  [1, 1, 1, 1, 1, 1, 1, 1]])
    Ainv = np.linalg.inv(A)
    p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n)
    p0s[-1] = 0
    p0s = np.asarray(p0s)
    deltas = np.array([delta(k, 0), delta(k, 1), delta(k, 2), delta(k, 3), delta(k, 4), delta(k, 5), delta(k, 6), 0])
    inhomgeneity = np.subtract(p0s, deltas)
    f_from_k = Ainv.dot(inhomgeneity)
    return f_from_k


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
        f = -(ca - 0.33) / tau + j * N * mean_puff_single(ca)
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


def intensities(N, j):
    Ds = []
    cas = np.linspace(0.01, 1.00, 100)
    for ca in cas:
        D = intensity_puff_single(ca)
        Ds.append((j ** 2) * N * D)
    return Ds


def derivation_intensities(N, j):
    dDs = []
    Ds = intensities(N, j)
    for i in range(100):
        if i == 0:
            dDs.append((Ds[i + 1] - Ds[i]) / 0.01)
        elif i == 99:
            dDs.append((Ds[i] - Ds[i - 1]) / 0.01)
        else:
            dDs.append((Ds[i + 1] - Ds[i - 1]) / 0.02)
    return dDs


def g(fs, dDs, Ds):
    gs = []
    for f, dD, D in zip(fs, dDs, Ds):
        g = (f - dD) / D
        gs.append(g)
    return gs


if __name__ == "__main__":
    home = os.path.expanduser("~")
    taua = 100
    ampa = 0.2
    n_cha = 4
    n_clu = 10

    taus = np.logspace(-1, 2, 30)
    js = np.logspace(-3, 0, 30)

    tau = taus[14]
    j = js[18]
    print(tau, j)
    data = np.loadtxt(
        home + "/CLionProjects/PhD/calcium_spikes_langevin/out/ca_langevin_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(
            tau, j, n_clu))
    ts, cas, jpuffs, adaps = np.transpose(data)

    # Calculate f(x) and D(x) from simulation. Convention dx/dt = f(x) + \sqrt{2D(x)}xi(t)
    # Hence f(x) = j_leak + <j_puff> and D(x) = j^2 * N * D_puff(x)
    jpuffs_over_ca = [[] for _ in range(100)]
    for ca, jpuff in zip(cas, jpuffs):
        for i in range(100):
            if ca > 0.01 * i and ca < 0.01 * (i + 1):
                jpuffs_over_ca[i].append(jpuff)
    mean_j_over_ca = []
    sqrt2D_j_over_ca = []
    cas_list = np.linspace(0, 1, 100)
    for ca, jpuff_over_ca in zip(cas_list, jpuffs_over_ca):
        mean_j_over_ca.append(np.mean(jpuff_over_ca) - (ca - 0.33) / tau)
        sqrt2D_j_over_ca.append(np.std(jpuff_over_ca))

    mean_add_std_sim = [mean + sqrt2D for mean, sqrt2D in zip(mean_j_over_ca, sqrt2D_j_over_ca)]
    mean_sub_std_sim = [mean - sqrt2D for mean, sqrt2D in zip(mean_j_over_ca, sqrt2D_j_over_ca)]

    # Calculate f(x) and D(x) from theory. Convention dx/dt = f(x) + \sqrt{2D(x)}xi(t)
    # Hence f(x) = j_leak + <j_puff> and D(x) = j^2 * N * D_puff(x)
    cas = np.linspace(0, 1, 100)
    fs = means(n_clu, tau, j)
    Ds = intensities(n_clu, j)  # Ds = j^2 * N * D_puff(x)
    dDs = derivation_intensities(n_clu, j)
    sqrt2D = [np.sqrt(2 * x) for x in Ds]
    drift_over_noises = g(fs, dDs, Ds)

    # ax.hist(cas, bins=100, density = True, alpha=.6, label="Langevin")
    # ax.hist(cas2, bins=100, density = True, alpha=.6, label="Markov")

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gs.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    remove_top_right_axis([ax])

    ax.plot(cas, mean_j_over_ca, c="k", label="f(x) sim")
    ax.plot(cas, fs, c="k", alpha=0.7, ls=":", label="f(x) theory")

    ax.plot(cas, sqrt2D_j_over_ca, label="$\sqrt{2D(x)}$ sim", c="C0", )
    ax.plot(cas, sqrt2D, label="$\sqrt{2D(x)}$ theory", c="C0", alpha=0.7, ls=":")

    ax.axvline(0.33, c="C7", ls="--")
    ax.plot(cas, dDs, label="$\partial_x D(x)$ theory", ls=":", c="C3")

    ax.legend(ncol=2)
    # ax.set_ylabel(r"$P_0(x)$")
    ax.set_xlabel(r"$[\rm{Ca}^{2+}]_i$")
    plt.savefig(home + "/Data/Calcium/Plots/8_langevin_drift_diffusion_coef.pdf", transparent=True)
    plt.show()
