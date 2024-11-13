import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st

if __name__ == "__main__":
    parameters = [[10.2, 8.82e-03, 4.14e+02, 1.78e-01],
                 [14.1, 9.90e-03, 6.49e+02, 5.78e-02],
                 [11.6, 1.16e-02, 1.94e+03, 4.53e-02],
                 [5.34, 2.14e-02, 5.12e+02, 1.44e-01],
                 [4.20, 1.86e-02, 6.50e+02, 3.68e-02],
                 [12.9, 1.19e-02, 8.02e+02, 1.60e-01],
                 [2.71, 2.53e-02, 1.75e+02, 8.52e-02],
                 [5.48, 1.47e-02, 7.59e+02, 6.30e-02],
                 [5.48, 1.47e-02, 5.01e+02, 9.45e-02],
                 [5.48, 1.47e-02, 9.81e+02, 1.90e-02],
                 [52.1, 6.80e-03, 6.00e+02, 1.15e-01],
                 [6.46, 1.51e-02, 5.79e+02, 5.23e-02],
                 [11.4, 1.13e-02, 9.36e+02, 5.84e-02],
                 [2.43, 3.13e-02, 2.48e+02, 6.63e-02],
                 [12.2, 8.48e-03, 1.54e+02, 8.71e-02],
                 [2.78, 2.67e-02, 3.10e+02, 6.25e-02],
                 [3.93, 2.05e-02, 1.60e+03, 4.10e-02],
                 [10.3, 8.65e-03, 4.15e+02, 5.27e-02],
                 [18.6, 6.02e-03, 9.31e+02, 4.84e-02],
                 [18.6, 6.02e-03, 9.40e+02, 4.83e-02],
                 [7.56, 1.17e-02, 3.65e+02, 8.35e-02],
                 [15.5, 7.47e-03, 1.92e+03, 3.41e-02],
                 [7.71, 1.09e-02, 5.37e+02, 2.38e-02],
                 [6.94, 1.74e-02, 1.19e+03, 5.05e-02]]

    tau_ers = [par[2] for par in parameters]
    eps_ers = [par[3] for par in parameters]

    st.set_default_plot_style()
    w  = 3.25*1.25
    h = 1.5*1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    axis = [ax1, ax2]
    st.remove_top_right_axis(axis)
    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')

    ax1.set_ylabel(r"$p(\tau_{\rm er})$")
    ax1.set_xlabel(r"$\tau_{\rm er}$")
    ax1.set_xlim([0, 2000])
    # ax1.set_xlim([100, 3000])
    # ax1.set_xscale("log")
    ax1.set_yticklabels([])
    ax1.hist(tau_ers, density=True, bins=[0, 400, 800, 1200, 1600, 2000], color=st.colors[4])
    ax1.axvline(np.mean(tau_ers), c="k", label=r"$\mu(\tau_{\rm er})$")
    ax1.axvline(np.mean(tau_ers) - np.std(tau_ers), c="k", ls=":", label=r"$\mu(\tau_{\rm er}) \pm \sigma(\tau_{\rm er})$")
    ax1.axvline(np.mean(tau_ers) + np.std(tau_ers), c="k", ls=":")
    # ax1.legend(fancybox=False, fontsize=8)

    ax2.set_ylabel(r"$p(\varepsilon)$")
    ax2.set_xlabel(r"$\varepsilon$")
    ax2.set_xlim([0, 0.2])
    # ax2.set_xlim([0.01, 1])
    # ax2.set_xscale("log")
    ax2.set_yticklabels([])
    ax2.hist(eps_ers, density=True, bins=[0, 0.04, 0.08, 0.12, 0.16, 0.2], color=st.colors[4])
    ax2.axvline(np.mean(eps_ers), c="k", label=r"$\mu(\varepsilon)$")
    ax2.axvline(np.mean(eps_ers) - np.std(eps_ers), c="k", ls=":", label=r"$\mu(\varepsilon) \pm \sigma(\varepsilon)$")
    ax2.axvline(np.mean(eps_ers) + np.std(eps_ers), c="k", ls=":")
    ax2.legend(fancybox=False, fontsize=8)

    plt.show()
