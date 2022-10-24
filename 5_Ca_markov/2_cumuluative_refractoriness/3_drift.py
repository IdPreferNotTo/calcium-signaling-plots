import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 3))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    axis = [ax1, ax2]
    st.remove_top_right_axis(axis)

    taus = [5, 1]
    ps = [0.015, 0.06]
    tau_er = 100
    eps_er = 0.1
    K = 10
    N = 5
    M = 3
    ip3 = 1.0
    cR = 0.2
    cT = 0.5
    for tau, p, ax in zip(taus, ps, axis):
        ax.set_xlabel(r"$c_{\rm i}$")
        ax.set_ylabel(r"$f(c_{\rm i}; c_{\rm er})$")
        ax.axhline(0, c="C7")
        adaps_fix = np.logspace(-2, 0, 5)
        cers_fix = 1 - adaps_fix
        for cer_fix in cers_fix:
            cis = np.linspace(cR, cT, 100)
            fs = []
            for ci in cis:
                f = fc.drift_theory(ci, p, tau, K, N, M, 1.0, cer=cer_fix)
                fs.append(f)
            ax.plot(cis, fs, label=f"{cer_fix:.2f}")
    legend = ax.legend(fancybox=False, framealpha=1., edgecolor="k", fontsize=8)
    plt.show()
