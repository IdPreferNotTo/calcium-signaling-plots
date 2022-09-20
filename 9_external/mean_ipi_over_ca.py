import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

import functions as fc
import styles as st

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 2.5))
    gs = gs.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    st.remove_top_right_axis([ax1, ax2])

    c0 = 0.2
    ropn_at_c0 = 0.5
    ropn_max = ropn_at_c0 * (1 + np.power(c0, 3))/np.power(c0, 3)
    rref = 20.0
    N = 5
    M = 3
    cas = np.linspace(0.2, 0.5, 100)
    means = []
    means_ref = []
    means_opn = []
    for ca in cas:
        ropn = ropn_max * (np.power(ca, 3) / (1. + np.power(ca, 3)))
        mean = (M - 1) / rref + 1 / ropn
        mean1 = (M - 1) / rref
        mean2 = 1/ropn
        means.append(mean)
        means_ref.append(mean1)
        means_opn.append(mean2)
    ax1.plot(cas, means, color=st.colors[0])
    ax1.plot(cas, means_ref, ls=":", color=st.colors[0])
    ax1.plot(cas, means_opn, ls="--", color=st.colors[0])
    ax1.set_xlabel(r"$c_{\rm i}$")
    ax1.set_ylabel(r"$\langle I \rangle$")
    ax1.set_xticks([0.2, 0.5])
    ax1.set_xticklabels(["$c_0$", "$c_T$"])

    means_dif_ref = []
    means_dif_ref1 = []
    means_dif_ref2 = []
    rrefs = [1, 10, 100]
    for rref in rrefs:
        means = []
        means_ref = []
        means_opn = []
        for ca in cas:
            ropn_at_c0 = 1/(2.1 - (M-1)/rref)
            ropn_max = ropn_at_c0 * (1 + np.power(c0, 3)) / np.power(c0, 3)
            ropn = ropn_max * (np.power(ca, 3) / (1. + np.power(ca, 3)))
            mean = (M - 1) / rref + 1 / ropn
            mean1 = (M - 1) / rref
            mean2 = 1 / ropn
            means.append(mean)
            means_ref.append(mean1)
            means_opn.append(mean2)
        means_dif_ref.append(means)
        means_dif_ref1.append(means_ref)
        means_dif_ref2.append(means_opn)
    for i, (means, means_ref, means_opn) in enumerate(zip(means_dif_ref, means_dif_ref1, means_dif_ref2)):
        ax2.plot(cas, means, color=st.colors[i+1])
        ax2.plot(cas, means_ref, ls=":", color=st.colors[i+1])
    plt.show()
