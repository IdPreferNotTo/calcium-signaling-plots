import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy import symbols, Eq, solve
from LiRinzel import utilities as ut

ip3 = 0.3
ca0, h0 = ut.fixed_point(ip3)
cas = []
dcas = []
dhs = []
ca = 0.05
h = np.cbrt(ut.h3_nullcline(ca, ip3))
for _ in range(10000):
    dca, dh = ut.f([ca, h], ip3)
    dcas.append(dca)
    dhs.append(dh)
    cas.append(ca)
    ca += dca*0.01
    h += dh*0.01

ca_crit, h_crit = ut.critical_point(ip3)
cas2 = []
dcas2 = []
dhs2 = []
ca = ca_crit -0.01
h = h_crit
for _ in range(1000):
    dca, dh = ut.f([ca, h], ip3)
    dcas2.append(dca)
    dhs2.append(dh)
    cas2.append(ca)
    ca += dca*0.01
    h += dh*0.01

f, ax = plt.subplots(1, 1, figsize=(6, 9 / 2))
ax.plot(cas, dcas, c="C0", label="$dca/dt$")
#ax.plot(cas2, dcas2)
#ax.plot(cas, dhs, c="C0", ls="--", label="$dh/dt$")

ax.axhline(0, ls="--", c="k", zorder=3)
ax.legend()
ax.set_ylabel(r"$d$Ca$^{2+}/dt$")
ax.set_xlabel(r"Ca$^{2+}$")
plt.show()

