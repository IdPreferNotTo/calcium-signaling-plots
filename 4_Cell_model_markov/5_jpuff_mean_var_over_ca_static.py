import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, TransformedBbox, BboxPatch, BboxConnector
import os


def coarse_grained_list(list, f):
    coarse_grained_list = []
    for i in range(0, len(list), f):
        average = sum(list[i:i+f])/f
        coarse_grained_list.append(average)
    return coarse_grained_list


def delta(a, b):
    if a==b:
        return 1
    else:
        return 0


def steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n):
    A = np.array([[-n*r_ref, 0, 0, 0, 0, 0, 0, r_cls],
                  [n*r_ref, -n*r_ref, 0, 0, 0, 0, 0, 0],
                  [0, n*r_ref, -n*r_ref, 0, 0, 0, 0, 0],
                  [0, 0, n*r_ref, -n*r_opn, 0, 0, 0, 0],
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
    deltas = np.array([delta(k, 0), delta(k, 1), delta(k, 2), delta(k, 3), delta(k,4), delta(k, 5), delta(k, 6), 0])
    inhomgeneity = np.subtract(p0s, deltas)
    f_from_k = Ainv.dot(inhomgeneity)
    return f_from_k


def gauss_dist(xs, mean, std):
    gauss_dist = []
    for x in xs:
        gauss = 1/np.sqrt(2*np.pi*(std**2)) * np.exp(-((x - mean)**2)/(2*std**2))
        gauss_dist.append(gauss)
    return gauss_dist


if __name__ == "__main__":
    #cas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    cas = np.linspace(0.01, 0.99, num=99)
    means_sim = []
    means_theory = []
    vars_sim = []
    vars_theory = []

    jpuffs_dt1_ca03 = []
    jpuffs_dt1_ca09 = []
    for ca in cas:
        print(ca)
        ca_fix = ca
        ca_rest = 0.33
        N = 10
        m = 4
        n = 4
        meanI = 10
        k=0
        home = os.path.expanduser("~")
        folder = "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/"
        file = "ca_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N{:d}_0.dat".format(ca_fix, N)
        data_spike = np.loadtxt(home + folder + file)

        ts, cals, jpuffs, adap = np.transpose(data_spike)
        means_sim.append(np.mean(jpuffs))
        dt0 = 0.1
        f = 10
        cg_list = coarse_grained_list(jpuffs, f)
        vars_sim.append(np.var(cg_list) * dt0 * f)

        r_opn = 0.13 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
        r_ref = 1.3 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
        r_cls = 50

        p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n)

        xs = [0, 0, 0, 0, 4, 3, 2, 1]
        idxs = [0, 1, 2, 3, 4, 5, 6, 7]
        mean = sum([N*x*p for x, p in zip(xs, p0s)])
        means_theory.append(mean)
        D_theory = 0
        for k in idxs:
            sum_over_i = 0
            f_from_k_to = f_from_k_invert_A(k, r_ref, r_opn, r_cls, m, n)
            for i in idxs:
                sum_over_i += xs[i]*f_from_k_to[i]
            D_theory += xs[k] * p0s[k]*sum_over_i
        vars_theory.append(2*N*D_theory)

        if ca == 0.3:
            jpuffs_dt1_ca03 = coarse_grained_list(jpuffs, 10)
        if ca == 0.9:
            jpuffs_dt1_ca09 = coarse_grained_list(jpuffs, 10)

    stds_sim = [np.sqrt(vars) for vars in vars_sim]
    stds_theory = [np.sqrt(vars) for vars in vars_theory]

    std_upper_theory = [mean + std for mean, std in zip(means_theory, stds_theory)]
    std_lower_theory = [mean - std for mean, std in zip(means_theory, stds_theory)]
    std_upper_sim= [mean + std for mean, std in zip(means_sim, stds_sim)]
    std_lower_sim = [mean - std for mean, std in zip(means_sim, stds_sim)]

    fig = plt.figure(tight_layout=True, figsize=(64/9, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    ax = fig.add_subplot(gs[0:2, 0:2])
    ax1 = fig.add_subplot(gs[0, 2])
    ax2 = fig.add_subplot(gs[1, 2])
    ax.plot(cas, means_sim, c="k", label =r"$\mu$")
    ax.plot(cas, means_theory, c="k", ls="--", label=r"$\mu$ theory")
    ax.fill_between(cas, std_upper_sim, std_lower_sim, color="C0", alpha=0.55, label=r"$\sigma_{y_1}\sqrt{\Delta t}, \Delta t = 1$")
    ax.plot(cas, std_upper_theory, c="C7", ls="--", label=r"$\sqrt{2D_N}$ theory")
    ax.plot(cas, std_lower_theory, c="C7", ls="--")
    ax.axhline(0, lw="1", c="k")
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    ax.set_xlabel(r"$[\rm{Ca}^{2+}]/K_{\rm Ca}$")
    ax.set_ylabel(r"$\mu, \sqrt{2D_N}$")
    ax.set_xlim([0, 1])
    ax.legend(loc=4, fancybox=False, framealpha=1.0)

    ax1.hist(jpuffs_dt1_ca09, bins=20, color="C0", alpha=0.55, density=True,
                  label=r"$c = 0.9$")
    ax1.set_xlabel("$y_1$")
    ax1.set_ylabel("$P(y_1)$")
    ax1.legend()
    mean_ca09 = np.mean(jpuffs_dt1_ca09)
    std_ca09 = np.std(jpuffs_dt1_ca09)
    rect = patches.Rectangle((0.89, mean_ca09 - 1.5*std_ca09), 0.02, 3*std_ca09, linewidth=1, edgecolor='k', facecolor='none')
    ax.add_patch(rect)

    gauss_ca09 = np.linspace(mean_ca09 - 3 * std_ca09, mean_ca09 + 3 * std_ca09, 100)
    ax1.plot(gauss_ca09, gauss_dist(gauss_ca09, mean_ca09, std_ca09), c="C7", ls="--")
    #ax2.legend(fancybox=False, framealpha=1.0, loc=4)

    ax2.hist(jpuffs_dt1_ca03, bins=20, color="C0", alpha=0.55, density=True,
                  label=r"$c = 0.3$")
    ax2.set_xlabel("$y_1$")
    ax2.set_ylabel("$P(y_1)$")
    ax2.legend()
    mean_ca03 = np.mean(jpuffs_dt1_ca03)
    std_ca03 = np.std(jpuffs_dt1_ca03)
    rect = patches.Rectangle((0.29, mean_ca03 - 1.5*std_ca03), 0.02, 3*std_ca03, linewidth=1, edgecolor='k', facecolor='none')
    ax.add_patch(rect)
    gauss_ca03 = np.linspace(mean_ca03 - 3 * std_ca03, mean_ca03 + 3 * std_ca03, 100)
    ax2.plot(gauss_ca03, gauss_dist(gauss_ca03, mean_ca03, std_ca03), c="C7", ls="--")

    axins = inset_axes(ax, width="30%", height="30%", loc="upper left", bbox_to_anchor=(0.2,0,1,1), bbox_transform=ax.transAxes)
    axins.plot(cas, stds_sim, c="C0", alpha=0.7)
    axins.plot(cas, stds_theory, c="C7", ls=":")
    axins.tick_params(direction="in")
    axins.set_xticks([0, 0.5, 1])
    #axins.set_yticks([0, 0.2, 0.4, 0.6])
    axins.set_xlabel(r"$[\rm{Ca}^{2+}]/K_{\rm Ca}$")
    axins.set_ylabel(r"$\sigma_{y_1}\sqrt{\Delta t}$")
    plt.savefig(home + "/Data/Calcium/Plots/markov_jpuff_mean_var_over_ca_static.pdf", transparent=True)
    plt.show()