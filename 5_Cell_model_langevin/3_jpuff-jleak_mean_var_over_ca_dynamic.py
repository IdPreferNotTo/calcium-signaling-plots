import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import functions as fc
import styles as st

if __name__ == "__main__":
    home = os.path.expanduser("~")
    taus = np.logspace(-1, 2, 50)
    jcas = np.logspace(-3, 0, 50)
    N = 10
    m = 4
    n = 5
    for i in range(1):
        tau = taus[30]
        j = jcas[25]
        print(tau, j)
        folder = "/CLionProjects/PhD/calcium_spikes_langevin/out/Data_no_adap/"
        file_ca = f"ca_langevin_ip1.00_tau{tau:.2e}_j{j:.2e}_N{N:d}_0.dat"
        file_spike = f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{j:.2e}_N{N:d}_0.dat"
        data = np.loadtxt(home + folder + file_ca)
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

        ISIs = np.loadtxt(home + folder + file_spike)
        mean_ISI = np.mean(ISIs)
        cv2_ISI = np.var(ISIs) / (mean_ISI ** 2)
        print(mean_ISI, cv2_ISI)

        st.set_default_plot_style()
        fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
        gs = gridspec.GridSpec(1, 1)
        ax = fig.add_subplot(gs[:])
        st.remove_top_right_axis([ax])

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
            meanI = 10
            r_opn = 0.13 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
            r_ref = 1.3 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
            r_cls = 50

            p0s = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)
            xs = fc.get_states(n, m)
            idxs = [i for i in range(n+m)]
            mean = sum([N*x*p for x, p in zip(xs, p0s)])
            means_theory.append(j*mean - (ca - 0.33)/tau)

            D_theory = 0
            for k in idxs:
                sum_over_i = 0
                f_from_k_to = fc.f_from_k_invert_M(k, r_ref, r_opn, r_cls, n, m)
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
        plt.savefig(home + "/Data/Calcium/Plots/langevin_puff_current_mean_var_dynamic.pdf", transparent=True)
        plt.show()
