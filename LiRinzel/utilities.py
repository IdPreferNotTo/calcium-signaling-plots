import numpy as np
from scipy.optimize import fsolve, fmin
import collections

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

def inverted_h3_nullcline(ca, ip3):
    m_inf = (ip3 / (ip3 + d1)) * (ca / (ca + d5))
    ca_er = (c0 - ca) / c1
    return -(-v2 / (v1 * np.power(m_inf, 3)) - v3 * np.power(ca, 2) / (
                (k3 ** 2 + ca ** 2) * (c1 * v1 * np.power(m_inf, 3) * (ca - ca_er))))


def h3_nullcline(ca, ip3):
    m_inf = (ip3 / (ip3 + d1)) * (ca / (ca + d5))
    ca_er = (c0 - ca) / c1
    return -v2 / (v1 * np.power(m_inf, 3)) - v3 * np.power(ca, 2) / (
                (k3 ** 2 + ca ** 2) * (c1 * v1 * np.power(m_inf, 3) * (ca - ca_er)))


def f(ca_h, ip3):
    Q2 = d2 * (ip3 + d1) / (ip3 + d3)
    ca, h = ca_h
    m_inf = (ip3 / (ip3 + d1)) * (ca / (ca + d5))
    h_inf = Q2 / (Q2 + ca)
    tau_h = 1. / (a2 * Q2 + a2 * ca)
    ca_er = (c0 - ca) / c1

    dca = -c1 * v1 * np.power(m_inf, 3) * np.power(h, 3) * (ca - ca_er) - c1 * v2 * (ca - ca_er) \
          - v3 * np.power(ca, 2) / (np.power(k3, 2) + np.power(ca, 2))
    dh = (h_inf - h) / tau_h
    return [dca, dh]


def fixed_point(ip3):
    ca_fix, h_fix = fsolve(f, [0.5, 0.5], args=ip3)
    if isinstance(ca_fix, collections.Sequence):
        ca_fix = min(ca_fix)
        h_fix = max(h_fix)
    return ca_fix, h_fix


def nullclines(ca0, h0, ip3):
    ca_null = np.linspace(0.01, ca0 + 0.5, 100)
    h_null = np.linspace(0.011, h0 + 0.5, 100)
    ca_nullcline = []
    h_nullcline = []
    for ca0, h0 in zip(ca_null, h_null):
        Q2 = d2 * (ip3 + d1) / (ip3 + d3)
        m_inf = (ip3 / (ip3 + d1)) * (ca0 / (ca0 + d5))
        h_inf = Q2 / (Q2 + ca0)
        tau_h = 1. / (a2 * Q2 + a2 * ca0)
        ca_er = (c0 - ca0) / c1
        h_nullcline.append((-v2 / (v1 * np.power(m_inf, 3)) - v3 * np.power(ca0, 2) / (
                (k3 ** 2 + ca0 ** 2) * (c1 * v1 * np.power(m_inf, 3) * (ca0 - ca_er)))) ** (1 / 3))
        ca_nullcline.append(Q2 / h0 - Q2)
    return [[ca_null, h_nullcline], [ca_nullcline, h_null]]


def time_series(ip3):
    t = 0
    dt = 0.0001
    ca = 0
    h = 0
    cat = []
    ht = []
    while t < 50:
        cat.append(ca)
        ht.append(h)
        ca_tmp = ca
        h_tmp = h
        dca, dh = f([ca_tmp, h_tmp], ip3)
        ca += dca * dt
        h += dh * dt
        t += dt
    return [cat, ht]


def max_nullcline_values(ip3):
    maximum = fmin(inverted_h3_nullcline, x0=0.1, args=(ip3,), full_output=True, disp=False)
    print(maximum)
    ca_max = maximum[0][0]
    h_max = -maximum[1]
    return ca_max, h_max ** (1 / 3)


def critical_point(ip3):
    # The passing of the hillock of the Calcium nullcline is a good criteria for the occurence of a spike in the excitabble
    # regime. In the tonically firing regime however it is not. Whats a spike in a oscillatory system anyway. However,
    # to find "spike times" in the oscillatory regime the criteria is ca > ca_fixed where ca_fixed is the ca coordinate
    # of the unstable fixed point
    ca_max, h_max = max_nullcline_values(ip3)
    ca_fix, h_fix = fixed_point(ip3)
    if (ca_max < ca_fix):
        ca_crit = ca_fix
        h_crit = h_fix
    else:
        ca_crit = ca_max
        h_crit = h_max
    return ca_crit, h_crit


def reset_values(ts, cas, h3s, ip3):
    ca_crit_start, h_crit_start = critical_point(ip3)
    h3_crit_start = h_crit_start ** 3
    start = []
    stop = []
    puffing = 0

    hs_reset = []
    cas_reset = []
    for t, c, h3 in zip(ts, cas, h3s):
        # if n_open > n_open_crit_start and c > ca_crit_start and puffing == 0:
        if c > ca_crit_start and puffing == 0:
            start.append(t)
            puffing = 1
        if c < ca_crit_start and puffing == 1 and h3 > h3_nullcline(c, ip3):
            stop.append(t)
            hs_reset.append(h3)
            cas_reset.append(c)
            puffing = 0

    return [cas_reset, hs_reset]


def spike_times(ts, cas, h3s, ip3):
    ca_crit_start, h_crit_start = critical_point(ip3)
    h3_crit_start = h_crit_start ** 3
    start = []
    stop = []
    puffing = 0

    for t, c, h3 in zip(ts, cas, h3s):
        #if n_open > n_open_crit_start and c > ca_crit_start and puffing == 0:
        if c > ca_crit_start and puffing == 0:
            start.append(t)
            puffing = 1
        if c < ca_crit_start and puffing == 1 and h3 > h3_nullcline(c, ip3):
            stop.append(t)
            puffing = 0

    ipis = []
    for t1, t2 in zip(start[:-1], start[1:]):
        ipis.append(t2 - t1)
    return ipis, start, stop


def released_ca(cas, n_opens, ip3):
    ca_crit, h_crit = critical_point(ip3)
    n_open_crit = h_crit ** 3
    areas = []
    area = 0
    puffing = 0
    for c, n_open in zip(cas, n_opens):
        if n_open > n_open_crit and c > ca_crit and puffing == 0:
            puffing = 1
        if n_open < n_open_crit and c < ca_crit and puffing == 1:
            puffing = 0
            areas.append(area)
            area = 0
        if puffing == 1:
            area += (c - ca_crit) * 0.05
    return areas


def k_corr(data1, data2, k):
    # Get two data set and calculate their correlation with lag k.
    mean1 = np.mean(data1)
    data1 = [x - mean1 for x in data1]
    mean2 = np.mean(data2)
    data2 = [x - mean2 for x in data2]
    k_corr = 0
    if k == 0:
        for x, y in zip(data1, data2):
            k_corr += x * y
    else:
        for x, y in zip(data1[:-k], data2[k:]):
            k_corr += x * y
    return k_corr / len(data2[k:])


def k_corr_var(data1, data2, k):
    # Get two arbitrary data set and calculate their correlation with lag k.
    k_corr = []
    if k == 0:
        for x, y in zip(data1, data2):
            k_corr.append(x * y)
    else:
        for x, y in zip(data1[:-k], data2[k:]):
            k_corr.append(x * y)
    k_corr_var = np.var(k_corr)
    return k_corr_var
