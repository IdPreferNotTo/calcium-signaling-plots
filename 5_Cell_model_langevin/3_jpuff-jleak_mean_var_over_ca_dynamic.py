import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
    taus = np.logspace(-1, 2, 50)
    jcas = np.logspace(-3, 0, 50)
    taua = 100
    ampa = 0.2
    n_cha = 4
    n_clu = 10
    for i in range(1):
        tau = taus[30] # [taus[30], taus[35], taus[40]][i]
        j = jcas[25] # [jcas[15], jcas[25], jcas[35]][i]
        print(tau, j)
        data = np.loadtxt(home + "/CLionProjects/PhD/calcium_spikes_langevin/out/ca_langevin_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu))
        ts, cas, jpuffs, adaps = np.transpose(data)

        jpuffs_over_ca = [[] for _ in range(100)]
        jleaks_over_ca = [[] for _ in range(100)]
        for ca, jpuff in zip(cas, jpuffs):
            for i in range(100):
                if ca > 0.01*i and ca < 0.01*(i+1):
                    jpuffs_over_ca[i].append(jpuff)

        mean_jpuff_over_ca = []
        std_jpuff_over_ca = []
        cas_list = np.linspace(0, 1, 100)
        for ca, jpuff_over_ca in zip(cas_list, jpuffs_over_ca):
            mean_jpuff_over_ca.append(np.mean(jpuff_over_ca) - (ca - 0.33)/tau)
            std_jpuff_over_ca.append(np.std(jpuff_over_ca))
        std_upper_sim = [mean + std for mean, std in zip(mean_jpuff_over_ca, std_jpuff_over_ca)]
        std_lower_sim = [mean - std for mean, std in zip(mean_jpuff_over_ca, std_jpuff_over_ca)]

        ISIs = np.loadtxt(
            home + "/CLionProjects/PhD/calcium_spikes_langevin/out/spike_times_langevin_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(
                tau, j, n_clu))

        mean_ISI = np.mean(ISIs)
        cv2_ISI = np.var(ISIs) / (mean_ISI ** 2)

        fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
        gs = gridspec.GridSpec(1, 1)
        ax = fig.add_subplot(gs[:])
        cas = np.linspace(0.0, 1, 100)
        ax.plot(cas, mean_jpuff_over_ca, label="$ISI$= {:.2f} \n CV = {:.2f}".format(mean_ISI, np.sqrt(cv2_ISI)))
        ax.fill_between(cas, std_upper_sim, std_lower_sim, color="C0", alpha=0.55, label=r"$\sigma_\tau, \tau = 1$")


        # Plot Theory
        means_theory = []
        stds_theory = []
        cas = np.linspace(0.01, 0.99, num=99)
        for ca in cas:
            ca_fix = ca
            ca_rest = 0.33
            N = 10
            m = 4
            n = 4
            meanI = 10
            r_opn = 0.13 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
            r_ref = 1.3 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
            r_cls = 50

            p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n)

            xs = [0, 0, 0, 0, 4, 3, 2, 1]
            idxs = [0, 1, 2, 3, 4, 5, 6, 7]
            mean = sum([N*x*p for x, p in zip(xs, p0s)])
            means_theory.append(j*mean - (ca - 0.33)/tau)

            D_theory = 0
            for k in idxs:
                sum_over_i = 0
                f_from_k_to = f_from_k_invert_A(k, r_ref, r_opn, r_cls, m, n)
                for i in idxs:
                    sum_over_i += xs[i] * f_from_k_to[i]
                D_theory += xs[k] * p0s[k] * sum_over_i
            stds_theory.append(j * np.sqrt(2 * N * D_theory))

        ax.plot(cas, means_theory, c="k")
        std_upper_theory = [mean + std for mean, std in zip(means_theory, stds_theory)]
        std_lower_theorry = [mean - std for mean, std in zip(means_theory, stds_theory)]
        ax.fill_between(cas, std_upper_theory, std_lower_theorry, color="C77", alpha=0.50, label=r"$\sqrt{2D}$")
        ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
        ax.set_xlabel(r"$[\rm{Ca}^{2+}]/K_{\rm Ca}$")
        ax.set_ylabel(r"$\mu, \sigma$")
        ax.set_xlim([0, 1])
        ax.legend()
        plt.savefig(home + "/Data/Calcium/Plots/langevin_variance_D_varying_ca_all_nClu1_nCha4_adap_rcls50.pdf",
                    transparent=True)
        plt.show()
