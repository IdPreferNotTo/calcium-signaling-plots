import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut

if __name__ == "__main__":
    ip3 = 0.3
    N =32

    data = np.loadtxt(
        "/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(N, ip3),
        skiprows=10_000)
    ts, cas, j1, h3s = np.transpose(data)
    cas_reset, hs_reset = ut.reset_values(ts, cas, h3s, ip3)

    f, ax = plt.subplots(1, 1, figsize=(6, 9 / 2))


    print("mean:", np.mean(cas_reset), np.mean(hs_reset))
    print("var:", np.var(cas_reset), np.var(hs_reset))
    plt.hist(cas_reset)
    plt.show()

    corr = []
    var = ut.k_corr(cas_reset, cas_reset, 0)
    ks = range(1,6)
    for k in ks:
        corr.append(ut.k_corr(cas_reset, cas_reset, k)/var)

    plt.plot(ks, corr)
    plt.show()
