import collections
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import fsolve, fmin
from LiRinzel import utilities as ut

for i in range(1,2):
    ip3 = 0.1 * (i + 1)
    ca_fix, h_fix = ut.fixed_point(ip3)

    [ca0, h_nullcline], [ca_nullcline, h0] = ut.nullclines(ca_fix, h_fix, ip3)

    N = 8
    data = np.loadtxt(
        "/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(N, ip3))
    t, ca, j1, n_open = np.transpose(data)

    ca_max, h_max = ut.max_nullcline_values(ip3)
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(ca0, h_nullcline, c="C3", ls="--")
    ax.plot(ca_nullcline, h0, c="C3", ls="--")
    ax.set_xlabel(r"$Ca{2+}$")
    ax.set_ylabel("$h$")
    ax.scatter(ca_fix, h_fix, c="k", zorder=3)
    ax.scatter(ca_max, h_max, c="k", zorder=3)
    ax.set_xlim(0, ca_fix + 0.5)
    ax.set_ylim(0, h_fix + 0.5)
    ca = ca[200:10000]
    n_open = n_open[200:10000]
    n_open = [x ** (1 / 3) for x in n_open]
    ax.plot(ca, n_open, c="C0")
    plt.savefig("/home/lukas/Plots/lr_ca_phd/stochastic_lr_phase_portrait_N{:d}_ip{:.2f}.pdf".format(N, ip3),
                transparent=True)
    plt.show()
    plt.close()
