import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut

def mean_ipi_var(ip3, N):
    data = np.loadtxt("/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(N, ip3), skiprows=10_000)
    ts, ca, j1, n_open = np.transpose(data)

    ipis, starts, stops = ut.spike_times(ts, ca, n_open, ip3)

    r0mean = np.mean([1/ipi for ipi in ipis])
    r0var = np.var([1/ipi for ipi in ipis])
    return r0mean, r0var

if __name__ == "__main__":

    ip3 = 0.3
    ip3s = [0.1*(i) for i in range(2, 10)]
    Ns = [2**(i) for i in range(1, 7)]
    N = 64
    r0means=[]
    r0vars=[]
    for ip3 in ip3s:
        r0mean, r0var = mean_ipi_var(ip3, N)
        r0means.append(r0mean)
        r0vars.append(r0var)

    f, ax = plt.subplots(1, 1, figsize=(6, 9 / 2))
    ax.errorbar(ip3s, r0means, yerr=r0vars)
    plt.show()