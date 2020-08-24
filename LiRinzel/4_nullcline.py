import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut


cas_null =  np.linspace(0., 1, 400)
hs_null = np.linspace(0., 1, 400)

ca0 = []
h0 = []
h1 = []
h2 = []

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

ip3 = 0.3
for h in hs_null:
    Q2 = d2 * (ip3 + d1) / (ip3 + d3)
    ca0.append(Q2 / h - Q2)

for ca in cas_null:
    h0.append(np.cbrt(ut.h3_nullcline(ca, ip3)))

caf, hf = ut.fixed_point(ip3)
cas = []
hs = []
dcas = []
dhs = []
ca = caf*1.1
h = hf*1.1
for _ in range(10000):
    dca, dh = ut.f([ca, h], ip3)
    dcas.append(dca)
    dhs.append(dh)
    cas.append(ca)
    hs.append(h)
    ca += dca*0.01
    h += dh*0.01

plt.plot(cas, hs, c="k")
plt.plot(cas_null, h0, c="C0")
plt.plot(ca0, hs_null, c="C1")
plt.plot()
plt.ylabel("h")
plt.xlabel("ca")
plt.ylim(0., 1.)
plt.xlim(0.0, 0.2)
plt.savefig("/home/lukas/Plots/lr_ca_phd/nullcine_h_ip{:.2f}.pdf".format(ip3), transparent=True)
plt.show()