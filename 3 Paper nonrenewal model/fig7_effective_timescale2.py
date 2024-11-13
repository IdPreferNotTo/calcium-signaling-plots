import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st
import functions as fc
import default_parameters as df


def exponential_cer(t, cer0, cer8, tau):
    return cer8 + (cer0 - cer8) * np.exp(-t / tau)

if __name__ == "__main__":
    # Set Plot style
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 1.75 * 1.25
    fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(w, h))
    ax1 = axs[0]
    ax2 = axs[1]

    ax1.text(0.1, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    st.remove_top_right_axis([ax1, ax2])
    ax1.set_ylabel(r"$\tau_{\rm eff}$ / s")
    ax1.set_xlabel(r"$\tau_{er}$ / s")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_ylabel(r"$\tau_{\rm eff}$ / s")
    ax2.set_xlabel(r"$\varepsilon$")
    ax2.set_xscale("log")
    ax2.set_ylim([0, 360])

    # Use mean-driven parameters only.
    tau = 5
    p = 0.015
    eps_er_fix = 0.03
    tau_er_fix = 300
    tau_ers = np.logspace(1, 3, 11)
    eps_ers = np.logspace(-2, 0, 11)
    tau_effs = []
    tau_1s = []
    tau_2s = []

    # Calculate tau1, tau2 (approx to tau_eff) as a function of tau_er
    for tau_er in tau_ers:
        data = df.load_adaptation_markov_transient(tau, p, tau_er, eps_er_fix)
        mean_cer = np.mean(data, axis=0)
        tau_eff = fc.measure_tau_eff(tau, p, tau_er, eps_er_fix)
        tau_effs.append(tau_eff)
        tau_1_1, tau_1_2 = fc.calculate_tau_1(tau, p, tau_er, eps_er_fix)
        tau_2_1, tau_2_2 = fc.calculate_tau_2(tau, p, tau_er, eps_er_fix)
        tau_1s.append(tau_1_1)
        tau_2s.append(tau_2_1)

    # Plot tau_eff, tau1, tau2 over tau_er
    ax1.set_xlim([8, 1200])
    ax1.axvspan(xmin=8, xmax=30, color="C7", alpha=0.2)
    ax1.scatter(tau_ers, tau_effs, fc="w", ec=st.colors[1], s=20, label=r"$\tau_{\rm eff}$")
    ax1.plot(np.linspace(10, 1000, 10_000), np.linspace(10, 1000, 10_000), ls=(0, (1, 1)), c="C7", label=r"$\tau_{\rm er}$")
    #ax1.plot(tau_ers, tau_1s, c=st.colors[7], ls=":", label=r"$\tau_1$")
    ax1.plot(tau_ers, tau_2s, c=st.colors[1], ls=(0, (3, 1)), label=r"$\tau_{\rm theo}$")

    # Plot legend to span both axes
    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.36, 0.0), loc=3,
                     ncol=4, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(1.)


    tau_effs = []
    tau_1s = []
    tau_2s = []
    # Calculate tau1, tau2 (approx to tau_eff) as a function of eps_er
    for eps_er in eps_ers:
        data = df.load_adaptation_markov_transient(tau, p, tau_er_fix, eps_er)
        mean_cer = np.mean(data, axis=0)
        tau_eff = fc.measure_tau_eff(tau, p, tau_er_fix, eps_er)
        tau_effs.append(tau_eff)
        tau_1_1, tau_1_2 = fc.calculate_tau_1(tau, p, tau_er_fix, eps_er)
        tau_2_1, tau_2_2 = fc.calculate_tau_2(tau, p, tau_er_fix, eps_er)
        tau_1s.append(tau_1_1)
        tau_2s.append(tau_2_1)

    # Plot tau_eff, tau1, tau2 over eps_er
    ax2.scatter(eps_ers, tau_effs, fc="w", ec=st.colors[1], s=20)
    #ax2.plot(eps_ers, tau_1s, c=st.colors[7], ls=":")
    ax2.plot(eps_ers, tau_2s, c=st.colors[1], ls=(0, (3, 1)))
    ax2.plot([0.01, 1], [tau_er_fix, tau_er_fix], ls=(0, (1, 1)), c="C7")

    # Save and show figure
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig7.pdf", dpi=300, transparent=True)
    plt.show()


