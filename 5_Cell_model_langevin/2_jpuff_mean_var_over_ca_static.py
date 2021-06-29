import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import os

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


if __name__ == "__main__":
    home = os.path.expanduser("~")
    cas_fix = np.linspace(0.01, 0.99, num=99)
    means = []
    stds = []
    means_t = []
    stds_t = []
    for ca in cas_fix:
        print(ca)
        data = np.loadtxt(
        home + "/CLionProjects/PhD/calcium_spikes_langevin/out/fixed calcium/ca_langevin_cafix{:.2f}_taua1.00e+02_ampa2.00e-01_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca))
        ts, cas, jpuffs, adaps = np.transpose(data)
        mean = np.mean(jpuffs)
        std = np.std(jpuffs)
        means.append(mean)
        stds.append(std)

        r_opn = 0.13 * np.power(ca / 0.33, 3) * ((1 + 0.33 ** 3) / (1 + ca ** 3))
        r_ref = 1.3 * np.power(ca / 0.33, 3) * ((1 + 0.33 ** 3) / (1 + ca ** 3))
        r_cls = 50

        p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, 4, 4)

        xs = [0, 0, 0, 0, 4, 3, 2, 1]
        idxs = [0, 1, 2, 3, 4, 5, 6, 7]
        mean = sum([10 * x * p for x, p in zip(xs, p0s)])
        means_t.append(mean)

        D_theory = 0
        for k in idxs:
            sum_over_i = 0
            f_from_k_to = f_from_k_invert_A(k, r_ref, r_opn, r_cls, 4, 4)
            for i in idxs:
                sum_over_i += xs[i]*f_from_k_to[i]
            D_theory += xs[k] * p0s[k]*sum_over_i
        stds_t.append(np.sqrt(2*10*D_theory))


    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gs.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])

    std_upper = [mean + std for mean, std in zip(means, stds)]
    std_lower = [mean - std for mean, std in zip(means, stds)]

    std_t_upper = [mean + std for mean, std in zip(means_t, stds_t)]
    std_t_lower = [mean - std for mean, std in zip(means_t, stds_t)]

    ax.plot(cas_fix, means)
    ax.plot(cas_fix, means_t, c="k")
    ax.fill_between(cas_fix, std_upper, std_lower, ls="--", color="C0", alpha=0.55, label=r"$\sigma_\tau, \tau = 1$")
    ax.plot(cas_fix, std_t_upper, c="C7", ls="--", label=r"$\sqrt{2D}$ theory")
    ax.plot(cas_fix, std_t_lower, c="C7", ls="--")

    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    ax.set_xlabel(r"$[\rm{Ca}^{2+}]/K_{\rm Ca}$")
    ax.set_ylabel(r"$\mu, \sigma$")
    ax.set_xlim([0, 1])
    ax.legend(loc=4, fancybox=False, framealpha=1.0)
    plt.savefig(home + "/Data/Calcium/Plots/langevin_variance_D_fixed_ca_all_nClu1_nCha4_adap_rcls50.pdf",
                transparent=True)
    plt.show()

