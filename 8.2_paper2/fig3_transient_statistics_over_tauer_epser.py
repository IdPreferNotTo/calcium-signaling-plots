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
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3.25*1.25, 2.5*1.25))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axis_top = [ax1, ax2]
    axis_bot = [ax3, ax4]
    axis = axis_top + axis_bot
    st.remove_top_right_axis(axis)

    ax1.set_ylabel(r"$n_{\rm tr}$")
    ax1.set_xlabel(r"$\tau_{er}$")
    ax1.set_xscale("log")

    ax2.set_ylabel(r"$\Delta T$")
    ax2.set_xlabel(r"$\tau_{er}$")
    ax2.set_xscale("log")

    ax3.set_ylabel(r"$n_{\rm tr}$")
    ax3.set_xlabel(r"$\varepsilon_{er}$")
    ax3.set_xscale("log")

    ax4.set_ylabel(r"$\Delta T$")
    ax4.set_xlabel(r"$\varepsilon_{er}$")
    ax4.set_xscale("log")
    # Parameters
    taus = [5.0, 1.0]
    ps = [0.015, 0.06]
    eps_er_fix = 0.1
    tau_er_fix = 100
    tau_ers = np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)

    home = os.path.expanduser("~")
    for tau, p, in zip(taus, ps):
        dTs = []
        n_trs = []
        n_trs_theory = []
        n_trs_theory2 = []
        dTs_theory = []
        for tau_er in tau_ers:
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er_fix)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=[100, 150, 2])
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs.append(T8 - T0)
            n_trs.append(n_tr)
            n_trs_theory.append(-1/np.log((1-eps_er_fix)*np.exp(-T8/tau_er)))
            eps_er_fix_h = eps_er_fix / (1 - eps_er_fix/2)
            n_trs_theory2.append((tau_er/T8)/(1 + eps_er_fix_h * tau_er / T8))


            if tau == 5.0:
                r0_no_adap = 2.92e-02
            else:
                r0_no_adap = 2.73e-02
            isi_no_adap = 1/r0_no_adap
            r0s_file = home + f"/Data/calcium_spikes_theory/r0_selfcon_tau{tau:.2e}_j{p:.2e}_{eps_er_fix:.2e}_{tau_er:.2e}.dat"
            data = np.loadtxt(r0s_file)
            r0, r0self = np.transpose(data)
            r0_intercept = interpolated_intercept(r0, r0, r0self)
            isi_self = 1/r0_intercept[0][0]
            dTs_theory.append(isi_self - isi_no_adap)

        ax1.scatter(tau_ers, n_trs, fc="w", ec=st.colors[0], alpha=0.75, s=20)
        ax2.scatter(tau_ers, dTs, fc="w", ec=st.colors[0], alpha=0.75, s=20)
        #ax1.plot(tau_ers, n_trs_theory, c=st.colors[5])
        ax2.plot(tau_ers, dTs_theory, lw=1, c=st.colors[5])
        print(n_trs_theory2)
        ax1.plot(tau_ers, n_trs_theory, lw=1, c=st.colors[5], zorder=3)

        dTs = []
        n_trs = []
        n_trs_theory = []
        n_trs_theory2 = []
        dTs_theory = []
        for eps_er in eps_ers:
            print(eps_er)
            data_isi = df.load_spike_times_markov_transient(tau, p, tau_er_fix, eps_er)
            rows, cols = data_isi.shape
            idx_max = cols
            idxs = np.arange(idx_max)
            means_Tidx = [np.mean(data_isi[:, idx]) for idx in idxs] # calculate the mean column wise
            popt, pcov = curve_fit(fc.exponential_Ti, idxs, means_Tidx, p0=(100, 150, 2))
            T0 = popt[0]
            T8 = popt[1]
            n_tr = popt[2]
            dTs.append(T8 - T0)
            n_trs.append(n_tr)
            n_trs_theory.append(-1 / (np.log(1 - eps_er) -T8 / tau_er_fix))
            eps_er_h = eps_er / (1 - eps_er/2)
            n_trs_theory2.append((tau_er_fix/T8)/(1 + eps_er_h * tau_er_fix / T8))

            if tau == 5.0:
                r0_no_adap = 2.92e-02
            else:
                r0_no_adap = 2.73e-02
            isi_no_adap = 1/r0_no_adap
            r0s_file = home + f"/Data/calcium_spikes_theory/r0_selfcon_tau{tau:.2e}_j{p:.2e}_{eps_er:.2e}_{tau_er_fix:.2e}.dat"
            data = np.loadtxt(r0s_file)
            r0, r0self = np.transpose(data)
            r0_intercept = interpolated_intercept(r0, r0, r0self)
            isi_self = 1/r0_intercept[0][0]
            dTs_theory.append(isi_self - isi_no_adap)

        ax3.scatter(eps_ers, n_trs, fc="w", ec=st.colors[0], alpha=0.75, s=20)
        ax4.scatter(eps_ers, dTs, fc="w", ec=st.colors[0], alpha=0.75, s=20)
        #ax3.plot(eps_ers, n_trs_theory, c=st.colors[5])
        ax4.plot(eps_ers, dTs_theory, lw=1, c=st.colors[5])
        ax4.set_ylim([0, 500])
        ax3.plot(eps_ers, n_trs_theory, lw=1, c=st.colors[5], zorder=3)

    home = os.path.expanduser("~")
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig3.pdf", transparent=True)
    plt.show()