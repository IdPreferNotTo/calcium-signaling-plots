import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ip3 = 0.30

def lr(ca, h):
    v1 = 6
    v2 = 0.11
    v3 = 0.9
    c0 = 2.
    c1 = 0.185
    k3 = 0.1
    a2 = 0.2
    d1 = 0.13
    d2 = 1.049
    d3 = 0.9434
    d5 = 0.08234

    Q2 = d2 * (ip3 + d1) / (ip3 + d3)
    m_inf = (ip3 / (ip3 + d1)) * (ca / (ca + d5))
    h_inf = Q2 / (Q2 + ca)
    tau_h = 1. / (a2 * Q2 + a2 * ca)
    ca_er = (c0 - ca) / c1

    dca = -c1 * v1 * np.power(m_inf, 3) * np.power(h, 3) * (ca - ca_er) - c1 * v2 * (ca - ca_er) \
          - v3 * np.power(ca, 2) / (np.power(k3, 2) + np.power(ca, 2))
    dh = (h_inf - h) / tau_h
    return dca, dh

ca_max, h_max = ut.max_nullcline_values(ip3)
ca = ca_max
h = h_max +0.02
cas = []
hs = []
ts = np.linspace(0, 20, 4000)
dt = ts[1]-ts[0]
for t in ts:
    cas.append(ca)
    hs.append(h)
    dca, dh = lr(ca, h)
    ca += dca*dt
    h += dh*dt

[ca0, h_nullcline], [ca_nullcline, h0] = ut.nullclines(0.5, 0.5, ip3)

fig, ax = plt.subplots(1,1)
ax.plot(ca0, h_nullcline, c="k")
ax.plot(ca_nullcline, h0, c="k")
ax.plot(cas, hs)
ax.set_ylim([min(hs)-0.05, max(hs)+0.05])
ax.set_xlim([min(cas)-0.05, max(cas)+0.05])
ax.set_ylabel(r"$h$")
ax.set_xlabel(r"Ca$^{2+}$")
ax_inset = inset_axes(ax, width=3, height=1, loc=1)
ax_inset.plot(ts, cas)
plt.show()