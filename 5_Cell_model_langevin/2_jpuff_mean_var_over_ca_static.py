import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import os

import functions as fc
import styles as st


if __name__ == "__main__":
    home = os.path.expanduser("~")
    cas_fix = np.linspace(0.01, 0.99, num=99)
    means = []
    stds = []
    means_t = []
    stds_t = []
    for ca in cas_fix:
        print(f"{ca:.2f}")
        data = np.loadtxt(
        home + "/CLionProjects/PhD/calcium_spikes_langevin/out/ca_fix/ca_langevin_cafix{:.2f}_ip1.00_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca))
        ts, cas, jpuff_mean, jpuff_std, adaps = np.transpose(data)
        mean = np.mean(jpuff_mean)
        std = np.std(jpuff_std)
        means.append(mean)
        stds.append(std)

        N = 10
        n = 5
        m= 4
        r_opn = 0.13 * np.power(ca / 0.33, 3) * ((1 + 0.33 ** 3) / (1 + ca ** 3))
        r_ref = 1.3 * np.power(ca / 0.33, 3) * ((1 + 0.33 ** 3) / (1 + ca ** 3))
        r_cls = 50

        p0s = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)

        xs = fc.get_states(n, m)
        idxs = [i for i in range(n+m)]
        mean = sum([N * x * p for x, p in zip(xs, p0s)])
        means_t.append(mean)

        D_theory = 0
        for k in idxs:
            sum_over_i = 0
            f_from_k_to = fc.f_from_k_invert_M(k, r_ref, r_opn, r_cls, n, m)
            for i in idxs:
                sum_over_i += xs[i]*f_from_k_to[i]
            D_theory += xs[k] * p0s[k]*sum_over_i
        stds_t.append(np.sqrt(2*N*D_theory))

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gs.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    st.remove_top_right_axis([ax])

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
    ax.set_xlabel(r"$c_i$")
    ax.set_ylabel(r"$\mu, \sigma$")
    ax.set_xlim([0, 1])
    ax.legend(loc=4, fancybox=False, framealpha=1.0)
    plt.savefig(home + "/Data/Calcium/Plots/langevin_puff_current_mean_var_static.pdf", transparent=True)
    plt.show()

