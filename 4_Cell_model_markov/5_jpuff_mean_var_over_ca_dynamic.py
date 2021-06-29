import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


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


def delta(a, b):
    if a==b:
        return 1
    else:
        return 0


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


def coarse_grained_list(list, f):
    coarse_grained_list = []
    for i in range(0, len(list), f):
        average = sum(list[i:i+f])/f
        coarse_grained_list.append(average)
    return coarse_grained_list


if __name__ == "__main__":
    home = os.path.expanduser("~")
    taus = np.logspace(-1, 2, 50)
    jcas = np.logspace(-3, 0, 50)
    Ncl = 10
    Nch = 4
    for i in range(1):
        tau = taus[30] # [taus[30], taus[35], taus[40]][i]
        jca = jcas[25] # [jcas[15], jcas[25], jcas[35]][i]
        print(tau, jca)
        data = np.loadtxt(home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_no_adap/spikes_tau{:.2e}_j{:.2e}_Ncha4_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(tau, jca))
        ts, cas, jpuffs, adaps = np.transpose(data)

        ts_clear = []
        cas_clear = []
        jpuffs_clear = []
        for t, ca, jpuff, adap in zip(ts, cas, jpuffs, adaps):
            if ca!=1:
                ts_clear.append(t)
                cas_clear.append(ca)
                jpuffs_clear.append(jpuff)

        coarse_grain_factor = 10
        ts_dt1 = ts_clear[::coarse_grain_factor]
        cas_dt1 = cas_clear[::coarse_grain_factor]
        jpuffs_dt1 = coarse_grained_list(jpuffs_clear, coarse_grain_factor)

        jpuffs_over_ca = [[] for i in range(100)]
        for ca, jpuff in zip(cas_clear, jpuffs_clear):
            for i in range(100):
                if ca > 0.01*i and ca < 0.01*(i+1):
                    jpuffs_over_ca[i].append(jpuff)


        mean_jpuff_over_ca = []
        std_jpuff_over_ca = []
        for jpuff_over_ca in jpuffs_over_ca:
            mean_jpuff_over_ca.append(np.mean(jpuff_over_ca))
            std_jpuff_over_ca.append(np.std(jpuff_over_ca))
        std_upper_sim = [mean + std/np.sqrt(coarse_grain_factor) for mean, std in zip(mean_jpuff_over_ca, std_jpuff_over_ca)]
        std_lower_sim = [mean - std/np.sqrt(coarse_grain_factor) for mean, std in zip(mean_jpuff_over_ca, std_jpuff_over_ca)]

        t_spike = 0
        ca_fix = 0.10
        ISI = []
        ts_abs = []
        t_tmp = 0
        for t, Ca in zip(ts, cas):
            if Ca == 1:
                ISI.append(t - t_tmp)
                t_tmp = t
        mean_ISI = np.mean(ISI)
        cv2_ISI = np.var(ISI) / (mean_ISI ** 2)

        fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
        gs = gridspec.GridSpec(1, 1)
        ax = fig.add_subplot(gs[:])
        cas = np.linspace(0.0, 1, 100)
        ax.plot(cas, mean_jpuff_over_ca)
        ax.fill_between(cas, std_upper_sim, std_lower_sim, color="C0", alpha=0.55, label=r"$\sigma_{\Delta t}, \Delta t = 1$")


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
            means_theory.append(jca*mean)

            D_theory = 0
            for k in idxs:
                sum_over_i = 0
                f_from_k_to = f_from_k_invert_A(k, r_ref, r_opn, r_cls, m, n)
                for i in idxs:
                    sum_over_i += xs[i] * f_from_k_to[i]
                D_theory += xs[k] * p0s[k] * sum_over_i
            stds_theory.append(jca * np.sqrt(2 * N * D_theory))

        ax.plot(cas, means_theory, c="k")
        std_upper_theory = [mean + std for mean, std in zip(means_theory, stds_theory)]
        std_lower_theorry = [mean - std for mean, std in zip(means_theory, stds_theory)]
        ax.fill_between(cas, std_upper_theory, std_lower_theorry, color="C7", alpha=0.50, label=r"$\sqrt{2D}$")
        ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
        ax.set_xlabel(r"$[\rm{Ca}^{2+}]/K_{\rm Ca}$")
        ax.set_ylabel(r"$\mu, \sqrt{2D_N}$")
        ax.set_xlim([0, 1])
        ax.legend()
        plt.savefig(home + "/Data/Calcium/Plots/markov_jpuff_mean_var_over_ca_dynamic.pdf", transparent=True)

        plt.show()
