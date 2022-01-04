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
    folder = "/CLionProjects/PhD/calcium_spikes_langevin/out/"
    file_ca = "ca_langevin_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu)
    file_spikes = "spike_times_langevin_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu)
    ca_data = np.loadtxt(home + folder + file_ca)
    isi_data = np.loadtxt(home + folder + file_spikes)
    ts, cas, jpuffs, adaps = np.transpose(ca_data)

    # f is the whole deterministic dynamics including leak and mean puff current
    cas_plot = np.linspace(0, 1, 101, endpoint=True)

    p0s_sim_ca = [[] for _ in range(100)]
    fs_sim_ca =  []
    sqrt2ds_sim_ca = []

    for ca, jpuff in zip(cas, jpuffs):
        for i, (ca_min, ca_max) in enumerate(zip(cas_plot[:-1], cas_plot[1:])):
            if ca_min < ca and ca < ca_max:
                p0s_sim_ca[i].append(jpuff)

    for ca, p0_sim_ca in zip(cas_plot, p0s_sim_ca):
        if len(p0_sim_ca) == 0:
            fs_sim_ca.append(np.nan)
            sqrt2ds_sim_ca.append(np.nan)
        else:
            fs_sim_ca.append(np.mean(p0_sim_ca) - (ca - 0.33) / tau)
            sqrt2ds_sim_ca.append(np.std(p0_sim_ca))

    # Calculate f(x) and D(x) from theory. Convention dx/dt = f(x) + \sqrt{2D(x)}xi(t)
    # Hence f(x) = j_leak + <j_puff> and D(x) = j^2 * N * D_puff(x)

    sqrt2D_theo_ca = []
    f_theo_ca = []
    for ca in cas_plot:
        f_theo_ca.append(f_func(ca, tau, j, n_clu))
        sqrt2D_theo_ca.append(np.sqrt(2*d_func(ca, j, n_clu)))


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

    #mean_T = 0
    #integral = 0
    #for ca in cas_theory[1:]:
    #    print(ca)
    #    h = h_func(ca, tau, j, n_clu)
    #    d = d_func(ca, j, n_clu)
    #    if ca < 0.33:
    #        integral += (np.exp(h)/d)*dca
    #    else:
    #        mean_T += np.exp(-h)*integral*dca
    #        integral += (np.exp(h)/d)*dca

    #print(1/mean_T)

    norm = np.sum(p0s_theo_ca) * dca
    r0theo = 1/norm
    r0sim = 1/np.mean(isi_data)
    cvsim = np.std(isi_data)/np.mean(isi_data)
    print(r0theo)
    print(r0sim)
    p0s_theo_ca = [x / norm for x in p0s_theo_ca]

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

    label_ca_sim = "Simulation \n" + rf"$r_{{0}}$ = {r0sim:.2f}" + "\n" + rf"$C_{{V}}$ = {cvsim:.2f}"
    label_ca_theo = "Theory \n" + rf"$r_{{0}}$ = {r0theo:.2f}" + "\n" + rf"$C_{{V}}$ = "
    ax.plot(cas_theory[1:], p0s_theo_ca[::-1], lw=1, c="k", label=label_ca_theo)
    ax.hist(cas, bins=50, density=True, alpha=.6, color="C0",
            label=label_ca_sim)
    legend = ax.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    ax.set_ylabel(r"$P_0([\rm{Ca}^{2+}]_i)$")
    ax.set_xlabel(r"$[\rm{Ca}^{2+}]_i$ [a.u.]")
    plt.savefig(home + "/Data/Calcium/Plots/9_langevin_ca_steady_state_probability.pdf", transparent=True)
    plt.show()
