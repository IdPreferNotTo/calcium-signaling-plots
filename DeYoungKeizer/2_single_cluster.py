import numpy as np
import matplotlib.pyplot as plt

v1 = 6.
v2 = 0.11
v3 = 0.9

c0 = 2.
c1 = 0.185
k3 = 0.1

a1 = 400  # binding of any IP3
a2 = 0.2  # binding of inhibitory Ca2+
a5 = 20  # binding of activation Ca2+

b1 = 52.  # removal of IP3 if inhibitory Ca2 + is not bound
b2 = 0.2098  # removal of inhibitory Ca2 + if IP3 is bound
b3 = 377.2  # removal of IP3 if inhibitory Ca2 + is bound
b4 = 0.0289  # removal of inhibitory Ca2 + if IP3 is not bound
b5 = 1.6468  # removal of activation Ca2 +

ip3 = 2.0
ca = 0.1

t = 0
dt = 0.0001
ts = []
puffs = []

open = [True, True, False]
channel = []
subunit = [False, False, False]
for i in range(3):
    channel.append(list(subunit))

while t < 10:
    N_open_subunits = 0
    for subunit in channel:
        ip3_is_bound = subunit[0]
        caa_is_bound = subunit[1]
        cai_is_bound = subunit[2]

        # binding probabilities
        if not ip3_is_bound:
            p_ip3 = a1 * ip3
        # unbinding probabilities
        else:
            if not cai_is_bound:
                p_ip3 = b1
            else:
                p_ip3 = b3

        # binding probabilities
        if not caa_is_bound:
            p_caa = a5 * ca
        # unbinding probabilities
        else:
            p_caa = b5

        # binding probabilities
        if not cai_is_bound:
            p_cai = a2 * ca
        # unbinding probabilities
        else:
            if not ip3_is_bound:
                p_cai = b4
            else:
                p_cai = b2

        p = p_ip3 + p_caa + p_cai

        xi1 = np.random.uniform(0, 1, 1)
        if xi1 < p * dt:
            xi2 = np.random.uniform(0, 1, 1)
            if xi2 < p_ip3 / p:
                subunit[0] = not subunit[0]
            elif xi2 < (p_caa + p_ip3) / p:
                subunit[1] = not subunit[1]
            else:
                subunit[2] = not subunit[2]
        if subunit == open:
            N_open_subunits += 1
    if N_open_subunits >= 3:
        puffs.append(1)
        ts.append(t)
    else:
        puffs.append(0)
        ts.append(t)
    t+=dt
plt.plot(ts, puffs)
plt.show()
