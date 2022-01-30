import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import functions as fc
import styles as st

if __name__ == "__main__":
    taus = np.logspace(-1, 2, 50)
    jcas = np.logspace(-3, 0, 50)
    tau = taus[30]
    jca = jcas[25]

    print(tau, jca)
    ca_res = 0.33
    N = 10
    m = 4
    n = 5
    home = os.path.expanduser("~")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/Data_no_adap/"
    file = f"ca_markov_ip1.00_tau{tau:.2e}_j{jca:.2e}_N10_0.dat"
    file_spike = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{jca:.2e}_N{N:d}_0.dat"
    data = np.loadtxt(home + folder + file)
    ts, cas, jpuffs, adaps = np.transpose(data)

    ISIs = np.loadtxt(home + folder + file_spike)
    mean_ISI = np.mean(ISIs)
    cv2_ISI = np.var(ISIs) / (mean_ISI ** 2)
    print(mean_ISI, cv2_ISI)

    ts_clear = []
    cas_clear = []
    jpuffs_clear = []
    for t, ca, jpuff, adap in zip(ts, cas, jpuffs, adaps):
        if ca != 1:
            ts_clear.append(t)
            cas_clear.append(ca)
            jpuffs_clear.append(jpuff)

    coarse_grain_factor = 10
    ts_dt1 = ts_clear[::coarse_grain_factor]
    cas_dt1 = cas_clear[::coarse_grain_factor]
    jpuffs_dt1 = fc.coarse_grain_list(jpuffs_clear, coarse_grain_factor)

    jpuffs_over_ca = [[] for i in range(100)]
    for ca, jpuff in zip(cas_clear, jpuffs_clear):
        for i in range(100):
            if ca > 0.01 * i and ca < 0.01 * (i + 1):
                jpuffs_over_ca[i].append(jpuff)

    mean_jpuff_over_ca = []
    std_jpuff_over_ca = []
    for jpuff_over_ca in jpuffs_over_ca:
        mean_jpuff_over_ca.append(np.mean(jpuff_over_ca))
        std_jpuff_over_ca.append(np.std(jpuff_over_ca))
    std_upper_sim = [mean + std / np.sqrt(coarse_grain_factor) for mean, std in
                     zip(mean_jpuff_over_ca, std_jpuff_over_ca)]
    std_lower_sim = [mean - std / np.sqrt(coarse_grain_factor) for mean, std in
                     zip(mean_jpuff_over_ca, std_jpuff_over_ca)]

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

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    st.remove_top_right_axis([ax])
    cas = np.linspace(0.0, 1, 100)
    ax.plot(cas, mean_jpuff_over_ca)
    ax.fill_between(cas, std_upper_sim, std_lower_sim, color="C0", alpha=0.55,
                    label=r"$\sigma_J \Delta t$")

    # Plot Theory
    means_theory = []
    stds_theory = []
    cas = np.linspace(0.01, 0.99, num=99)
    for ca in cas:
        ca_fix = ca
        ca_rest = 0.33
        meanI = 10
        r_opn = 0.13 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest ** 3) / (1 + ca_fix ** 3))
        r_ref = 1.3 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest ** 3) / (1 + ca_fix ** 3))
        r_cls = 50

        p0s = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)
        xs = fc.get_states(n, m)

        mean = sum([N * x * p for x, p in zip(xs, p0s)])
        means_theory.append(jca * mean)

        D_theory = 0
        for k in range(n+m):
            sum_over_i = 0
            f_from_k_to =fc.f_from_k_invert_M(k, r_ref, r_opn, r_cls, n, m)
            for i in range(n+m):
                sum_over_i += xs[i] * f_from_k_to[i]
            D_theory += xs[k] * p0s[k] * sum_over_i
        stds_theory.append(jca * np.sqrt(2 * N * D_theory))

    ax.plot(cas, means_theory, c="k")
    std_upper_theory = [mean + std for mean, std in zip(means_theory, stds_theory)]
    std_lower_theorry = [mean - std for mean, std in zip(means_theory, stds_theory)]
    ax.fill_between(cas, std_upper_theory, std_lower_theorry, color="C7", alpha=0.50, label=r"$\sqrt{2D}$")
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    ax.set_xlabel(r"$c_i$")
    ax.set_ylabel(r"$\mu, \sqrt{2D_N}$")
    ax.set_xlim([0, 1])
    ax.legend()
    plt.savefig(home + "/Data/Calcium/Plots/puff_current_mean_var_dynamic.pdf", transparent=True)

    plt.show()
