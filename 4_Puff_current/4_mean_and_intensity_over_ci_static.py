import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    # Prepare data
    home = os.path.expanduser("~")
    means_sim = []
    means_theory = []

    vars_sim = []
    intensities_theory = []

    jps_dt1_cR = []
    jps_dt1_cT = []

    tau = 1.
    p = 1.
    cT = 0.50
    cR = 0.20
    ip3 = 1.
    N = 5
    M = 3
    ci = 0.20
    r_ref = df.r_ref
    r_cls = df.r_cls

    cis = np.logspace(-1, 0, 10)
    cis = [ci / 2 for ci in cis]
    for idx, ci in enumerate(cis):
        r_opn = df.r_opn(ci)
        print(f"{ci:.2f}")
        data = df.load_traces_fixed_ci_markov(tau, p, ci, K=1)
        ts, cas, jps, cers = np.transpose(data)
        dt = ts[1] - ts[0]
        f = int(1/dt)

        means_sim.append(np.mean(jps))
        cg_list = fc.coarse_grain_list(jps, f)
        vars_sim.append(np.var(cg_list) * dt * f) # dt * f = the actual time bin

        mean_theory = fc.mean_jp_single_theory(ci, N, M, s=ip3)
        intensity_theory = fc.noise_intensity_jp_single_theory(ci, N, M, s=ip3)

        means_theory.append(mean_theory)
        intensities_theory.append(intensity_theory)

        if idx == 5:
            jps_dt1_cR = fc.coarse_grain_list(jps, f)
        if idx == 9:
            jps_dt1_cT = fc.coarse_grain_list(jps, f)

    stds_sim = [np.sqrt(vars) for vars in vars_sim]
    stds_theory = [np.sqrt(2*D) for D in intensities_theory]

    std_upper_theory = [mean + std for mean, std in zip(means_theory, stds_theory)]
    std_lower_theory = [mean - std for mean, std in zip(means_theory, stds_theory)]
    std_upper_sim= [mean + std for mean, std in zip(means_sim, stds_sim)]
    std_lower_sim = [mean - std for mean, std in zip(means_sim, stds_sim)]

    # Plot data
    st.set_default_plot_style()
    fig = plt.figure(tight_layout = True, figsize=(64/9, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    ax = fig.add_subplot(gs[0:2, 0:2])
    ax1 = fig.add_subplot(gs[0, 2])
    ax2 = fig.add_subplot(gs[1, 2])
    st.remove_top_right_axis([ax, ax1, ax2])

    # Main ax
    ax.set_xlabel(r"$c_i$")
    ax.set_xlim([0, 0.5])
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_ylabel(r"$\mu, \sqrt{2D_N}$")
    ax.plot(cis, means_theory, c="k", ls="--", label=r"$\mu$ theory")
    ax.plot(cis, means_sim, c="k", label =r"$\mu$")
    ax.fill_between(cis, std_upper_sim, std_lower_sim, color="C0", alpha=0.55, label=r"$\sigma_{y_1}\sqrt{\Delta t}, \Delta t = 1$")
    ax.plot(cis, std_upper_theory, c="C0", ls="--", label=r"$\sqrt{2D_N}$ theory")
    ax.plot(cis, std_lower_theory, c="C0", ls="--")
    ax.axhline(0, lw="1", ls=":", c="k")
    legend = ax.legend(loc=2, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)

    # Histogram for ci = cR
    ax1.set_xlabel("$Y$")
    ax1.set_ylabel("$P(Y)$")
    ax1.hist(jps_dt1_cR, bins=50, color="C0", alpha=0.55, density=True)
    mean_cR = np.mean(jps_dt1_cR)
    std_cR = np.std(jps_dt1_cR)
    gauss_ca09 = np.linspace(mean_cR - 3 * std_cR, mean_cR + 3 * std_cR, 100)
    ax1.plot(gauss_ca09, fc.gaussian_dist(gauss_ca09, mean_cR, std_cR), c="C7")

    # Histogram for ca = cT
    ax2.set_xlabel("Y")
    ax2.set_ylabel("$P(Y)$")
    ax2.hist(jps_dt1_cT, bins=50, color="C0", alpha=0.55, density=True)
    mean_cT = np.mean(jps_dt1_cT)
    std_cT = np.std(jps_dt1_cT)
    gauss_ca03 = np.linspace(mean_cT - 3 * std_cT, mean_cT + 3 * std_cT, 100)
    ax2.plot(gauss_ca03, fc.gaussian_dist(gauss_ca03, mean_cT, std_cT), c="C7")

    # Save and show figure
    plt.show()