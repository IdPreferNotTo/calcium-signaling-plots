import collections
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import fsolve, fmin
from LiRinzel import utilities as ut

for i in range(10):
    ip3 = 0.1*(i+1)
    ca_fix, h_fix = fsolve(ut.f, [0.5, 0.5], args=ip3)
    if isinstance(ca_fix, collections.Sequence):
        ca_fix = min(ca_fix)
        h_fix = max(h_fix)
    ca = np.linspace(0.01, ca_fix+0.5, 20)
    h = np.linspace(0.01, h_fix+0.5, 20)
    CA, H = np.meshgrid(ca, h)
    u, v = np.zeros(CA.shape), np.zeros(H.shape)
    NI, NJ = CA.shape
    for i in range(NI):
        for j in range(NJ):
            ca = CA[i, j]
            h = H[i, j]
            yprime = ut.f([ca, h], ip3)
            u[i, j] = yprime[0]
            v[i, j] = yprime[1]

    [ca0, h_nullcline] , [ca_nullcline, h0]= ut.nullclines(ca_fix, h_fix, ip3)
    cat, ht = ut.time_series(ip3)

    fig = plt.figure()
    gs  = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    Q = ax.quiver(CA, H, u, v, color='k')
    ax.plot(cat[100000:], ht[100000:], c="C0")
    ax.plot(ca0, h_nullcline, c="C3", ls="--")
    ax.plot(ca_nullcline, h0, c="C3", ls="--")
    ax.set_xlabel(r"$Ca{2+}$")
    ax.set_ylabel("$h$")
    ax.set_xlim(0, ca_fix+0.5)
    ax.set_ylim(0, h_fix+0.5)
    plt.savefig("/home/lukas/Plots/lr_ca_phd/lr_phase_portrait_ip{:.2f}.pdf".format(ip3), transparent=True)
    plt.show()
    plt.close()