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


def exponential_cer(t, cer8, tau):
    return cer8 + (1 - cer8) * np.exp(-t / tau)

if __name__ == "__main__":
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.66 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1])

    ax1.text(0.05, 0.95, r"A$_1$", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.05, 0.95, r"A$_2$", fontsize=11, transform=ax2.transAxes, va='top')
    ax3.text(0.05, 0.95, r"B$_1$", fontsize=11, transform=ax3.transAxes, va='top')
    ax4.text(0.05, 0.95, r"B$_2$", fontsize=11, transform=ax4.transAxes, va='top')

    axis_top = [ax1, ax3]
    axis_bot = [ax2, ax4]
    axis = axis_top + axis_bot
    st.remove_top_right_axis(axis)


    ax1.set_ylabel(r"$n_{\rm tr}$")
    ax1.set_xlabel(r"$\tau_{er}$")
    ax1.set_xscale("log")

    ax2.set_ylabel(r"$\Delta T$")
    ax2.set_xlabel(r"$\tau_{er}$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax3.set_ylabel(r"$n_{\rm tr}$")
    ax3.set_xlabel(r"$\varepsilon$")
    ax3.set_xscale("log")

    ax4.set_ylabel(r"$\Delta T$")
    ax4.set_xlabel(r"$\varepsilon$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")

    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.05
    tau_er_fix = 500
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)

    home = os.path.expanduser("~")
    for i, (tau, p), in enumerate(zip(taus, ps)):
        dTs = []
        n_trs = []
        dTs_theory = []
        n_trs_theory = []
        n_trs_theory2 = []
        for tau_er in tau_ers:
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er_fix)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            func = lambda x, t8, ntr: fc.exponential_Ti(x, means_Tidx[0], t8, ntr)
            popt, pcov = curve_fit(func, idxs, means_Tidx, p0=(150, 2))
            T8 = popt[0]
            n_tr = popt[1]
            dTs.append(T8 - means_Tidx[0])
            n_trs.append(n_tr)
            eps_fix = eps_er_fix / (1 - eps_er_fix/2)
            n_trs_theory2.append((tau_er/T8)/(1 + eps_fix * tau_er / T8))

            T0_theory = fc.calculate_T_init(tau, p)
            T8_theory = fc.calculate_T_infty(tau, p, tau_er, eps_er_fix)
            dTs_theory.append(T8_theory - T0_theory)
            tau_1, tau_2 = fc.calculate_tau_1(tau, p, tau_er, eps_er_fix)
            n_trs_theory.append(tau_1/T0_theory)

        if i == 0:
            label="mean-driven"
            ax1.scatter(tau_ers, n_trs, fc="w", ec=st.colors[0], s=15, label=label)
            ax2.scatter(tau_ers, dTs, fc="w", ec=st.colors[0], s=15)
        else:
            label="excitable"
            ax1.scatter(tau_ers, n_trs, fc="w", ec=st.colors[1], s=15, label=label)
            ax2.scatter(tau_ers, dTs, fc="w", ec=st.colors[1], s=15)

        if i == 0:
            ax1.plot(tau_ers, n_trs_theory, c=st.colors[5], zorder=3, label="Theory")
        #else:
        #    ax1.plot(tau_ers, n_trs_theory, c=st.colors[5], zorder=3)
        ax1.plot(tau_ers, n_trs_theory2, c=st.colors[5], ls="--")
        ax2.plot(tau_ers, dTs_theory, c=st.colors[5])

        dTs = []
        n_trs = []
        dTs_theory = []
        n_trs_theory = []
        n_trs_theory2 = []
        for eps_er in eps_ers:
            print(eps_er)
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er_fix, eps_er)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            func = lambda x, t8, ntr: fc.exponential_Ti(x, means_Tidx[0], t8, ntr)
            popt, pcov = curve_fit(func, idxs, means_Tidx, p0=(150, 2))
            T8 = popt[0]
            n_tr = popt[1]
            dTs.append(T8 - means_Tidx[0])
            n_trs.append(n_tr)
            eps = eps_er / (1 - eps_er/2)

            n_trs_theory2.append((tau_er_fix/T8)/(1 + eps * tau_er_fix / T8))
            T0_theory = fc.calculate_T_init(tau, p)
            T8_theory = fc.calculate_T_infty(tau, p, tau_er_fix, eps_er)
            dTs_theory.append(T8_theory - T0_theory)
            tau_1, tau_2 = fc.calculate_tau_1(tau, p, tau_er_fix, eps_er)
            n_trs_theory.append(tau_1/T0_theory)

        if i == 0:
            ax3.scatter(eps_ers, n_trs, fc="w", ec=st.colors[0], s=15)
            ax4.scatter(eps_ers, dTs, fc="w", ec=st.colors[0], s=15)
            ax3.plot(eps_ers, n_trs_theory, c=st.colors[5])
        else:
            ax3.scatter(eps_ers, n_trs, fc="w", ec=st.colors[1], s=15)
            ax4.scatter(eps_ers, dTs, fc="w", ec=st.colors[1], s=15)
        ax3.plot(eps_ers, n_trs_theory2, c=st.colors[5], ls="--")

        ax4.plot(eps_ers, dTs_theory, c=st.colors[5])


    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.465, .0), loc=3,
                      ncol=4, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(1.)

    home = os.path.expanduser("~")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig7.pdf", transparent=True)
    plt.show()