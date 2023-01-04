from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

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
    tau = 5
    j = 0.015
    eps_er = 0.1
    tau_er = 100 #np.logspace(1, 3, 21)
    eps_ers = np.logspace(-2, 0, 21)
    for eps_er in eps_ers:
        home = os.path.expanduser("~")
        st.set_default_plot_style()
        fig = plt.figure(tight_layout=True, figsize=(8, 4))
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])
        st.remove_top_right_axis([ax1])

        ax1.set_xlabel("$r_0$")
        ax1.set_ylabel("$r_0$")

        r0s_file = home + f"/Data/calcium_spikes_theory/r0_selfcon_tau{tau:.2e}_j{j:.2e}_{eps_er:.2e}_{tau_er:.2e}.dat"
        data = np.loadtxt(r0s_file)
        r0, r0self = np.transpose(data)
        ax1.plot(r0, r0)
        ax1.plot(r0, r0self)

        data_isi = df.load_spike_times_markov(tau, j, cer=True, taua=tau_er, ampa=eps_er)
        mean_isi = np.mean(data_isi[20:])
        ax1.scatter(1/mean_isi, 1/mean_isi)

        x_inter, y_inter = interpolated_intercept(r0, r0, r0self)
        ax1.scatter(x_inter, y_inter)

        print(1/mean_isi, x_inter)
        plt.show()

