import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import fmin
from LiRinzel import utilities as ut



for j in range(10):
    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    axis = [ax1, ax2, ax3]
    for i in range(3):

        n_total = 2**(i+4)
        ip3 = 0.1*(j+1)
        print(n_total, ip3)
        data = np.loadtxt("/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(n_total, ip3))
        t, ca, j1, n_open = np.transpose(data)
        t_max = 200
        maximum = fmin(ut.h_nullcline, x0=0.1, args=(ip3,), full_output=True)
        print(maximum)
        ca_crit = maximum[0][0]
        h_crit = -maximum[1]
        axis[i].axhline(ca_crit, c="C0", ls="--")
        axis[i].axhline(h_crit, c="C2", ls="--")
        axis[i].set_xlim(400,600)
        axis[i].set_xlim(400,600)
        axis[i].plot(t, ca, label=r"$N_{{Channel}}=$ {:d}".format(n_total))
        axis[i].plot(t, n_open, lw=1, alpha=0.7, c="C2")
        axis[i].legend()
    plt.savefig("/home/lukas/Plots/lr_ca_phd/lr_timeseries_ip{:.2f}.pdf".format(ip3), transparent=True)
    plt.show()
    plt.close()