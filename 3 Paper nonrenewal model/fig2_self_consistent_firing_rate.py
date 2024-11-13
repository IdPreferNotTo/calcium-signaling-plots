import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st
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
    w = 3.25 * 1.25
    h = 1.75 * 1.25
    fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(w, h))
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.text(0.1, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')
    st.remove_top_right_axis([ax1, ax2])

    for i, (ax, tau, p) in enumerate(zip([ax1, ax2], [5, 1], [0.015, 0.06])):
        home = os.path.expanduser("~")

        # Load data to calculate firing rate r0 for fixed cer
        cers_m = np.linspace(0.7, 1, 100)
        data_fpe = np.loadtxt(home + f"/Data/calcium/theory/cer_r0_fpe_tau{tau:.2e}_j{p:.2e}.dat")

        cers_f, r0s_f = np.transpose(data_fpe)

        if i == 0:
            ax.plot(cers_f, r0s_f, color=st.colors[1], label=r"$r_0(\tau, p, \langle c_{\rm er} \rangle)$")
            ax.plot([1.1, 1.2], [1., 1.], color=st.colors[3], label=r"$r_0(\tau, p, \langle c_{\rm er} \rangle)$")
        else:
            ax.plot(cers_f, r0s_f, color=st.colors[3], label=r"$r_0(\tau, p, \langle c_{\rm er} \rangle)$")

        eps_er = 0.03
        tau_er = 300
        eps = eps_er/(1.-eps_er/2.)
        r0s_er = np.zeros(401)
        for i, cer in enumerate(cers_f):
            r0s_er[i] = (1. - cer) / (eps*tau_er*cer)
        ax.plot(cers_f, r0s_er, c="k", label=r"$r_0(\tau_{\rm er}, \varepsilon, \langle c_{\rm er} \rangle)$")

        x, y = interpolated_intercept(cers_f, r0s_er, r0s_f)
        ax.axhline(y, ls=":", c="C7")
        ax.axvline(x, ls=":", c="C7")

    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.5, .0), loc=3,
                     ncol=3, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(1.)
    ax1.set_ylim([-0.01, 0.05])
    ax1.set_xlim([0.8, 1])
    ax1.set_xlabel(r"$\langle c_{\rm er} \rangle$")
    ax1.set_ylabel(r"$r_0(\langle c_{\rm er} \rangle)$ / s\textsuperscript{-1}")

    ax2.set_ylim([-0.01, 0.05])
    ax2.set_xlim([0.8, 1])
    ax2.set_xlabel(r"$\langle c_{\rm er} \rangle$")
    ax2.set_ylabel(r"$r_0(\langle c_{\rm er} \rangle)$ /  s\textsuperscript{-1}")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig2.pdf", dpi=300, transparent=True)
    plt.show()