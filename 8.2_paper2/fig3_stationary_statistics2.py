import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

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


if __name__ == "__main__":
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25*1.25
    h = 4.0*1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(3, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[0, 1])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[2, 1])
    axis_top = [ax1, ax2, ax3]
    axis_bot = [ax4, ax5, ax6]
    axis = axis_top + axis_bot
    st.remove_top_right_axis(axis)
    ax1.text(0.15, 0.95, r"A$_1$", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.15, 0.95, r"A$_2$", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.15, 0.95, r"A$_3$", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.15, 0.95, r"B$_1$", fontsize=11, transform=ax4.transAxes, va='top')
    ax5.text(0.15, 0.95, r"B$_2$", fontsize=11, transform=ax5.transAxes, va='top')
    ax6.text(0.15, 0.95, r"B$_3$", fontsize=11, transform=ax6.transAxes, va='top')

    ax1.set_ylabel(r"$\langle T \rangle$")
    ax1.set_xlabel(r"$\tau_{er}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax2.set_ylabel(r"$C_V$")
    ax2.set_xlabel(r"$\tau_{er}$")
    ax2.set_xscale("log")
    ax2.set_ylim([0, 1])

    ax3.set_ylabel(r"$\rho_1$")
    ax3.set_xlabel(r"$\tau_{er}$")
    ax3.set_xscale("log")
    ax3.set_ylim([-0.5, 0.2])
    ax3.axhline(0, ls=":", c="C7")

    ax4.set_ylabel(r"$\langle T \rangle$")
    ax4.set_xlabel(r"$\varepsilon$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")

    ax5.set_ylabel(r"$C_V$")
    ax5.set_xlabel(r"$\varepsilon$")
    ax5.set_xscale("log")
    ax5.set_ylim([0, 1])

    ax6.set_ylabel(r"$\rho_1$")
    ax6.set_xlabel(r"$\varepsilon$")
    ax6.set_xscale("log")
    ax6.set_ylim([-0.5, 0.2])
    ax6.axhline(0, ls=":", c="C7")

    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.05
    tau_er_fix = 500
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)

    for i, (tau, p) in enumerate(zip(taus, ps)):
        means = []
        cvs = []
        p1s = []
        mean_isis_theory = []
        for tau_er in tau_ers:
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er_fix)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            cvs.append(cv)
            p1s.append(p1)

            mean_isi_theory = fc.calculate_T_infty(tau, p, tau_er, eps_er_fix)
            mean_isis_theory.append(mean_isi_theory)

        if i==0:
            label="mean-driven"
        else:
            label="excitable"
        ax1.scatter(tau_ers, means, fc="w", ec=st.colors[i], s=15, zorder=3, label=label)
        if i == 1:
            ax1.plot(tau_ers, mean_isis_theory, c=st.colors[5], zorder=3, label="Theory")
        else:
            ax1.plot(tau_ers, mean_isis_theory, c=st.colors[5], zorder=3)

        ax1.plot(tau_ers, tau_ers, c="C7", ls=":")
        ax2.scatter(tau_ers, cvs, fc="w", ec=st.colors[i], s=15, zorder=3)
        ax3.scatter(tau_ers, p1s, fc="w", ec=st.colors[i], s=15, zorder=3)

        means = []
        cvs = []
        p1s = []
        mean_isis_theory = []
        for eps_er in eps_ers:
            data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er_fix, ampa=eps_er)
            mean = np.mean(data_isi)
            std = np.std(data_isi)
            var = np.var(data_isi)
            cv = std/mean
            p1 = fc.k_corr(data_isi, data_isi, 1)/var
            means.append(mean)
            cvs.append(cv)
            p1s.append(p1)

            mean_isi_theory = fc.calculate_T_infty(tau, p, tau_er_fix, eps_er)
            mean_isis_theory.append(mean_isi_theory)

        ax4.scatter(eps_ers, means, fc="w", ec=st.colors[i], s=15, zorder=3)
        ax4.plot(eps_ers, mean_isis_theory, c=st.colors[5], zorder=3)
        ax5.scatter(eps_ers, cvs, fc="w", ec=st.colors[i], s=15, zorder=3)
        ax6.scatter(eps_ers, p1s, fc="w", ec=st.colors[i], s=15, zorder=3)

    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.465, .0), loc=3,
                      ncol=4, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(.75)
    home = os.path.expanduser("~")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig3.pdf", transparent=True)
    plt.show()