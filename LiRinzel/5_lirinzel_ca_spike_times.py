import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import fmin
from LiRinzel import utilities as ut


for j in range(10):
    fig = plt.figure()
    gs = gridspec.GridSpec(6, 1)
    ax2 = fig.add_subplot(gs[1:3, 0])
    ax1 = fig.add_subplot(gs[0, 0], sharex=ax2)
    ax4 = fig.add_subplot(gs[4:6, 0])
    ax3 = fig.add_subplot(gs[3, 0], sharex=ax4)
    axis = [ax1, ax2, ax3, ax4]
    for i in range(2):
        n_total = 10**(i+1)
        ip3 = 0.1*(j+1)
        print(n_total, ip3)
        data = np.loadtxt("/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(n_total, ip3))
        t, ca, j1, n_open = np.transpose(data)
        t_max = 200
        ca_crit, h_crit = ut.critical_value(ip3)
        n_open_crit = h_crit**3
        ipis, puff_times, puff_stop_times = ut.spike_times(t, ca, n_open, ip3)
        for puff_time, puff_stop_time in zip(puff_times, puff_stop_times):
            if puff_time > 400 and puff_time < 600:
                axis[2*i].axvline(puff_time)
        axis[2*i +1].axhline(ca_crit, c="C0", ls="--")
        axis[2*i +1].axhline(n_open_crit, c="C2", ls="--")
        axis[2*i +1].set_xlim(400,600)
        axis[2*i +1].set_xlim(400,600)
        axis[2*i +1].plot(t, ca, label=r"$N_{{Channel}}=$ {:d}".format(n_total))
        axis[2*i +1].plot(t, n_open, lw=1, alpha=0.7, c="C2")
        axis[2*i +1 ].legend()
    plt.savefig("/home/lukas/Plots/lr_ca_phd/lr_puff_times_ip{:.2f}.pdf".format(ip3), transparent=True)
    plt.show()
    plt.close()