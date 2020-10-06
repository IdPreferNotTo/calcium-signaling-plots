import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut


for i in range(1):
    ip3 = 0.3

    N = 32
    data = np.loadtxt("/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(N, ip3), skiprows=10_000)
    ts, ca, j1, n_open = np.transpose(data)

    ipis, starts, stops = ut.spike_times(ts, ca, n_open, ip3)

    plt.hist(ipis, bins=30, label="$IP_3 = {:.1f}$".format(ip3) +  "\n" +  r"$\langle \Delta t \rangle = {:.1f}$".format(np.mean(ipis)) + "\n" + "$\sigma = {:.1f}$".format(np.std(ipis)))
    plt.legend()
    plt.show()
