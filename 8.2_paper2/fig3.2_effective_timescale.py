import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st
import functions as fc


def interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    return xc,yc


def exponential_cer(t, cer0, cer8, tau):
    return cer8 + (cer0 - cer8) * np.exp(-t / tau)

if __name__ == "__main__":
    # Set Plot style
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 1.75 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax1.text(0.1, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    st.remove_top_right_axis([ax1, ax2])
    ax1.set_ylabel(r"$\tau_{\rm eff}$")
    ax1.set_xlabel(r"$\tau_{er}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_ylabel(r"$\tau_{\rm eff}$")
    ax2.set_xlabel(r"$\varepsilon$")
    ax2.set_xscale("log")
    ax2.set_ylim([0, 250])

    # Use mean-driven parameters only.
    tau = 5
    p = 0.015
    eps_er_fix = 0.05
    tau_er_fix = 200
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)
    tau_effs = []
    tau_0s = []
    tau_1s = []

    # Calculate tau1, tau2 (approx to tau_eff) as a function of tau_er
    for tau_er in tau_ers:
        file = f"transient_adaptation_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er_fix:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
        data = np.loadtxt("/home/lukas/CLionProjects/PhD/calcium/calcium_spikes_markov_transient/out/" + file,
                          usecols=np.arange(0, 1000))
        mean_cer = np.mean(data, axis=0)
        tau_eff = fc.measure_tau_eff(tau, p, tau_er, eps_er_fix)
        tau_effs.append(tau_eff)
        tau_1_1, tau_1_2 = fc.calculate_tau_0(tau, p, tau_er, eps_er_fix)
        tau_2_1, tau_2_2 = fc.calculate_tau_1(tau, p, tau_er, eps_er_fix)
        tau_0s.append(tau_1_1)
        tau_1s.append(tau_2_1)

    # Plot tau_eff, tau1, tau2 over tau_er
    ax1.set_xlim([8, 1200])
    ax1.axvspan(xmin=8, xmax=30, color="C7", alpha=0.5)
    ax1.scatter(tau_ers, tau_effs, fc="w", ec=st.colors[0], s=15, label=r"$\tau_{\rm eff}$")
    ax1.plot(np.linspace(10, 1000, 10_000), np.linspace(10, 1000, 10_000), ls=(0, (1, 1)), c="C7", label=r"$\tau_{\rm er}$")
    ax1.plot(tau_ers, tau_0s, c=st.colors[5], ls="--", label=r"$\tau_1$")
    ax1.plot(tau_ers, tau_1s, c=st.colors[5], label=r"$\tau_2$")

    # Plot legend to span both axes
    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.465, .0), loc=3,
                     ncol=4, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(1.)


    tau_effs = []
    tau_0s = []
    tau_1s = []
    # Calculate tau1, tau2 (approx to tau_eff) as a function of eps_er
    for eps_er in eps_ers:
        file = f"transient_adaptation_markov_ip{1.00:.2f}_taua{tau_er_fix:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K{10}_{5:d}.dat"
        data = np.loadtxt("/home/lukas/CLionProjects/PhD/calcium/calcium_spikes_markov_transient/out/" + file,
                          usecols=np.arange(0, 1000))
        mean_cer = np.mean(data, axis=0)
        tau_eff = fc.measure_tau_eff(tau, p, tau_er_fix, eps_er)
        tau_effs.append(tau_eff)
        tau_1_1, tau_1_2 = fc.calculate_tau_0(tau, p, tau_er_fix, eps_er)
        tau_2_1, tau_2_2 = fc.calculate_tau_1(tau, p, tau_er_fix, eps_er)
        tau_0s.append(tau_1_1)
        tau_1s.append(tau_2_1)

    # Plot tau_eff, tau1, tau2 over eps_er
    ax2.scatter(eps_ers, tau_effs, fc="w", ec=st.colors[0], s=15)
    ax2.plot(eps_ers, tau_0s, c=st.colors[5], ls="--")
    ax2.plot(eps_ers, tau_1s, c=st.colors[5])
    ax2.plot([0.01, 1], [tau_er_fix, tau_er_fix], ls=(0, (1, 1)), c="C7")

    # Save and show figure
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig3.2.pdf")
    plt.show()


