import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

for j in range(10):
    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    axis = [ax1, ax2, ax3]
    for i in range(3):
        n_total = 10**(i+1)
        ip3 = 0.1*(j+1)
        data = np.loadtxt("/home/lukas/CLionProjects/dyk-calcium-phd/out/dyk_model_ca_waves_N{:d}_ip{:.2f}".format(n_total, ip3))
        t, ca, j1, n_open = np.transpose(data)
        t_max = 200
        axis[i].set_xlim(500,600)
        axis[i].set_xlim(500,600)
        axis[i].plot(t, ca, label=r"$N_{{Channel}}=$ {:d}".format(n_total))
        axis[i].plot(t, n_open, lw=1, alpha=0.7, c="C2")
        axis[i].legend()
    plt.savefig("/home/lukas/Plots/dyk_ca_phd/dyk_ip{:.2f}.pdf".format(ip3), transparent=True)
    plt.show()
    plt.close()

