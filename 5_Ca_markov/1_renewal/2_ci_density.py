import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    home = os.path.expanduser("~")
    # Parameters
    tau = 5.
    p = 0.015
    ip3 = 1.0
    M = 3
    N = 5
    K = 10

    # Load Data
    data_ci = df.load_traces_markov(tau, p, cer=False)
    ts, cas, jpuffs, adaps = np.transpose(data_ci)

    data_theory = df.load_ci_density(tau, p)
    ci_theory, density_theory = np.transpose(data_theory)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    st.remove_top_right_axis([ax])

    ax.plot(ci_theory, density_theory, lw=1, c="k")
    ax.hist(cas, bins=50, density=True, alpha=.6, color="C0")
    ax.set_ylabel(r"$p_0(c_{\rm i})$")
    ax.set_xlabel(r"$c_{\rm i}$")
    plt.show()
