import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import os

import functions as fc
import styles as st


if __name__ == "__main__":
    home = os.path.expanduser("~")
    cas_fix = np.logspace(-1, 0, 10)
    means = []
    means_theo = []

    for ca in cas_fix:
        print(f"{ca:.2f}")
        data = np.loadtxt(
        home + "/Data/calcium_spikes_markov/ca_fix/ca_markov_cafix{:.2f}_ip1.00_tau1.00e+00_j1.00e+00_K1_5.dat".format(ca))
        ts, cas, jpuff_mean, adaps = np.transpose(data)
        mean = np.mean(jpuff_mean)
        means.append(mean)

        K = 1
        N = 5
        M = 3
        r_opn = 0.2
        r_ref = 20
        r_cls = 50
        fc.mean_puff_single(ca, N, M, 1, r_opn,  r_ref)

        means_theo.append(mean)



    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gs.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    st.remove_top_right_axis([ax])
    ax.scatter(cas_fix, means)
    ax.plot(cas_fix, means_theo)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    ax.set_xlabel(r"$c_i$")
    ax.set_ylabel(r"$\mu$")
    plt.savefig(home + "/Plots/puff_model_reduction_two_states_intensity.pdf", transparent=True)
    plt.show()

