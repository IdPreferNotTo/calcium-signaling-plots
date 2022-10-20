import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    tau = 5.
    p = 0.015
    ip3 = 1.0
    cR = 0.2
    K = 10
    M = 3
    N = 5
    data = df.load_traces_markov(tau, p, cer = False)
    ts, cas, jpuffs, adaps = np.transpose(data)
    dt = ts[1] - ts[0]

    ts_clear = []
    cas_clear = []
    jpuffs_clear = []
    for t, ca, jpuff, adap in zip(ts, cas, jpuffs, adaps):
        if ca != 0.5:
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
                jpuffs_over_ca[i].append(jpuff/p)

    mean_jpuff_over_ca = []
    std_jpuff_over_ca = []
    for jpuff_over_ca in jpuffs_over_ca:
        mean_jpuff_over_ca.append(np.mean(jpuff_over_ca))
        std_jpuff_over_ca.append(np.std(jpuff_over_ca))
    std_upper_sim = [mean + std / np.sqrt(coarse_grain_factor) for mean, std in
                     zip(mean_jpuff_over_ca, std_jpuff_over_ca)]
    std_lower_sim = [mean - std / np.sqrt(coarse_grain_factor) for mean, std in
                     zip(mean_jpuff_over_ca, std_jpuff_over_ca)]


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
    intensities_theory = []
    cis = np.linspace(0.01, 0.50, num=50)
    for ci in cis:
        mean_theory = fc.mean_jp_single_theory(ci, N, M, s=ip3)
        intensity_theory = fc.noise_intensity_jp_single_theory(ci, N, M, s=ip3)
        means_theory.append(mean_theory)
        intensities_theory.append(intensity_theory)

    stds_theory = [np.sqrt(2*K*D) for D in intensities_theory]
    means_theory = [K*mu for mu in means_theory]

    ax.plot(cis, means_theory, c="k")
    std_upper_theory = [mean + std for mean, std in zip(means_theory, stds_theory)]
    std_lower_theorry = [mean - std for mean, std in zip(means_theory, stds_theory)]
    ax.fill_between(cis, std_upper_theory, std_lower_theorry, color="C7", alpha=0.50, label=r"$\sqrt{2D}$")
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xlabel(r"$c_i$")
    ax.set_ylabel(r"$K\mu_x, \sqrt{2KD_x}$")
    ax.legend()
    plt.show()
