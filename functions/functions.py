import numpy as np
from typing import List
from scipy.integrate import quad

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def delta(a, b):
    """
    Classical delta function \delta(a-b)
    """
    if a == b:
        return 1
    else:
        return 0


def theta(a, b):
    if a < b:
        return 0
    else:
        return 1

def equally_spaced_puff_data(data, dt):
    # data = time, state, idx
    data = [x if x[1] > 0 else np.asarray([x[0], 0, x[2]]) for x in data]
    ts = np.arange(0, 4990, step=dt)
    data_equal = []
    n = 0
    for t in ts:
        if data[n][0] < t:
            n += 1
        if n == 0:
            data_equal.append([t, 0])
        else:
            data_equal.append([t, data[n-1][1]])
    return data_equal

def autocorrelation(data):
    mean = np.mean(data)
    data = [x - mean for x in data]
    corr = 0
    for x, y in zip(data[:-1], data[1:]):
        corr += x * y
    corr /= len(data[1:])
    return corr

def autocorrelation_function(data, dt):
    mean = np.mean(data)
    data = [x - mean for x in data]
    correlation = []
    ks = np.arange(0, 100)
    for k in ks:
        corr = 0
        if k == 0:
            for x, y in zip(data, data):
                corr += x * y
            corr /= len(data)
        else:
            for x, y in zip(data[:-k], data[k:]):
                corr += x * y
            corr /= len(data[k:])
        correlation.append(corr)
    return [dt*k for k in ks], correlation


def fano_factor_interspike_intervals(isis, dt):
    spike_time = 0
    spike_times = []
    for isi in isis:
        spike_time += isi
        spike_times.append(spike_time)
    maxt = sum(isis)

    bins = int(maxt/dt) + 1
    counts = np.zeros(bins)
    for spike_time in spike_times:
        counts[int(spike_time/dt)] += 1
    mean_count = np.mean(counts)
    var_count = np.var(counts)
    return mean_count, var_count


def puff_open_probability(data):
    data = [ele for ele in data if ele[2]==0]
    t, state, idx, sums = np.transpose(data)
    t_open = 0
    t_closed = 0
    for t1, t2, s1, s2 in zip(t[:-1], t[1:], state[:-1], state[1:]):
        if s1 <= 0:
            t_closed += t2 - t1
        else:
            t_open += t2 - t1
    return t_open/(t_open + t_closed)


def puff_tau_closed(data):
    data = [ele for ele in data if ele[2]==0]
    ts, ss, idx, sums = np.transpose(data)
    nr = 0
    for s in ss:
        if s == 1:
            nr += 1
    t_open = 0
    t_closed = 0
    for t1, t2, s1, s2 in zip(ts[:-1], ts[1:], ss[:-1], ss[1:]):
        if s1 <= 0:
            t_closed += t2 - t1
        else:
            t_open += t2 - t1
    return t_closed/nr

def puff_tau_open(data):
    data = [ele for ele in data if ele[2]==0]
    ts, ss, idx, sums = np.transpose(data)
    nr = -1
    for s in ss:
        if s == 1:
            nr += 1
    t_open = 0
    t_closed = 0
    for t1, t2, s1, s2 in zip(ts[:-1], ts[1:], ss[:-1], ss[1:]):
        if s1 <= 0:
            t_closed += t2 - t1
        else:
            t_open += t2 - t1
    return t_open/nr

def moments(xs, k):
    """
    Calculates the k-th moment of the sample data xs:
    1/n \sum^n (x - <x>)**k
    where n the the sample length.
    """
    moment = 0
    mu = np.mean(xs)
    for x in xs:
        moment += (x - mu)**k
    return moment/len(xs)


def gaussian_dist(xs, mean, std):
    ys = []
    for x in xs:
        y = 1/np.sqrt(2*np.pi*(std**2)) * np.exp(-((x - mean)**2)/(2*std**2))
        ys.append(y)
    return ys


def coarse_grain_list(l: List[float], f: int):
    """
    Create a coarse grained version of the original list where the elements of the new list
    are the mean of the previous list over some window that is determined by f.
    f determines how many elements are averaged l_new[i] = mean(l[f*(i):f*(i+1)])
    """
    l_new = []
    max = int(len(l)/f) - 1
    for i in range(max):
        mean = np.mean(l[f*i:f*(i+1)])
        l_new.append(mean)
    return l_new

def moving_coarse_grain_list(l: List[float], n: int):
    """
    Create a coarse grained version of the original list where the elements of the new list
    are the mean of the previous list over some window that is determined by f.
    f determines how many elements are averaged l_new[i] = mean(l[f*(i):f*(i+1)])
    """
    if n<1:
        return l
    else:
        l_new = []
        for i in range(n, len(l)-n):
            mean = np.mean(l[i-n:i+n])
            l_new.append(mean)
        return l_new


def get_states(n, m):
    xs = np.empty(n + m)
    for i in range(n + m):
        if i < m:
            xs[i] = 0
        if i >= m:
            xs[i] = n + m - i
    return xs


def steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m):
    # M is the actual transition matrix
    M = np.zeros([n+m, n+m])
    for i in range(n+m):
            if i < m-1:
                M[i, i] = -r_ref
                M[i+1, i] = r_ref
            elif i == m-1:
                M[i, i] = -r_opn
                for k in range(i+1, n+m):
                    M[k, i] = r_opn/n
            else:
                M[i, i] = -r_cls
                M[(i+1)%(n+m), i] = r_cls
    # M does not have full rank. Wherefor some row has to be replaced to reflect the additional normalization
    M[-1] = np.ones(n+m)
    Minv = np.linalg.inv(M)
    inhomgeneity = np.zeros(n+m)
    inhomgeneity[-1] = 1
    p0s = Minv.dot(inhomgeneity)
    return p0s


def f_from_k_invert_M(k, r_ref, r_opn, r_cls, n, m):
    M = np.zeros([n + m, n + m])
    for i in range(n + m):
        if i < m - 1:
            M[i, i] = - r_ref
            M[i + 1, i] = r_ref
        elif i == m - 1:
            M[i, i] = -r_opn
            for j in range(i + 1, n + m):
                M[j, i] = r_opn/n
        else:
            M[i, i] = -r_cls
            M[(i + 1) % (n + m), i] = r_cls
    M[-1] = np.ones(n + m)

    Minv = np.linalg.inv(M)
    p0s = steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)
    p0s[-1] = 0
    p0s = np.asarray(p0s)
    deltas = [delta(k, i) for i in range(n+m)]
    deltas[-1] = 0
    inhomgeneity = np.subtract(p0s, deltas)
    f_from_k = Minv.dot(inhomgeneity)
    return f_from_k


def mean_puff_single(x, n, m, IP3, r_opn_single, r_ref):
    r_opn = n * r_opn_single * np.power(x / 0.2, 3) * ((1 + 0.2 ** 3) / (1 + x ** 3)) * np.power(IP3 / 1., 3) * ((1. + 1. ** 3) / (1. + IP3 ** 3))
    r_ref = r_ref
    r_cls = 50

    p0s = steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)
    xs = get_states(n, m)
    mean = sum([x * p for x, p in zip(xs, p0s)])
    return mean


def intensity_puff_single(x, n, m, IP3, r_opn_single, r_ref):
    r_opn = n * r_opn_single * np.power(x / 0.20, 3) * ((1. + 0.20 ** 3) / (1. + x ** 3)) * np.power(IP3 / 1., 3) * ((1. + 1. ** 3) / (1. + IP3 ** 3))
    r_ref = r_ref #* np.power(x / 0.20, 3) * ((1. + 0.20 ** 3) / (1. + x ** 3)) * np.power(IP3 / 1., 3) * ((1. + 1. ** 3) / (1. + IP3 ** 3))
    r_cls = 50

    xs = get_states(n, m)
    idxs = [i for i in range(n+m)]
    p0s = steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)

    D_theory = 0
    for k in idxs:
        sum_over_i = 0
        f_from_k_to = f_from_k_invert_M(k, r_ref, r_opn, r_cls, n, m)
        for i in idxs:
            sum_over_i += xs[i] * f_from_k_to[i]
        D_theory += xs[k] * p0s[k] * sum_over_i
    return D_theory


def d_func(x, j, N, n, m, IP3, r_opn_single, r_ref):
    if x == 0:
        return 0
    else:
        return np.power(j, 2) * N * intensity_puff_single(x, n, m, IP3, r_opn_single, r_ref)


def f_func(x, tau, j, N, n, m, IP3, r_opn_single, r_ref):
    if x == 0:
        return -(x - 0.20) / tau
    else:
        f = mean_puff_single(x, n, m, IP3, r_opn_single, r_ref)
        return -(x - 0.20) / tau + j * N * f


def g_func(x, tau, j, N, n, m, IP3):
    f = f_func(x, tau, j, N, n, m, IP3)
    d = d_func(x, j, N, n, m, IP3)
    return f/d


def h_func(x, tau, j, N, n, m, IP3):
    #dca = 0.0001
    #h = 0
    #ca = 0.20
    #while(ca <= x):
    #    print(ca)
    #    g = g_func(ca, tau, j, N, n, m, IP3)
    #    h += g*dca
    #    ca += dca
    h = quad(g_func, 0.20, x, args=(tau, j, N, n, m, IP3))[0]
    return h


def firing_rate_no_adap(tau, j, N, n, m, IP3 = 1):
    cas_theory = np.linspace(0.30, 1, 10*700 + 1)
    dca = cas_theory[1] - cas_theory[0]
    p0s_theo_ca = []
    integral = 0

    for ca in reversed(cas_theory[1:]):
        print(f"{ca:.3f}")
        h = h_func(ca, tau, j, N, n, m, IP3)
        d = d_func(ca, j, N, n, m, IP3)
        if ca == 1:
            integral += 0
        elif ca >= 0.20:
            integral += np.exp(-h)*dca
        p0s_theo_ca.append(integral * np.exp(h) / d)
    print(p0s_theo_ca)
    norm = np.sum(p0s_theo_ca) * dca
    r0 = 1 / norm
    return r0

def k_corr(data1, data2, k):
    # Get two arbitrary data set and calculate their correlation with lag k.
    mu1 = np.mean(data1)
    mu2 = np.mean(data2)
    data1 = [x - mu1 for x in data1]
    data2 = [x - mu2 for x in data2]
    k_corr = 0
    if k == 0:
        for x, y in zip(data1, data2):
            k_corr += x * y
    else:
        for x, y in zip(data1[:-k], data2[k:]):
            k_corr += x * y
    return k_corr / len(data2[k:])


def fourier_transformation_isis(f, isis):
    t = 0
    trafo = 0
    for isi in isis:
        t += isi
        trafo += np.exp(1j*2*np.pi*f*t)
    return trafo

def power_spectrum_isis(fs, isis, Tmax=2000):
    ISIs_chunks = []
    chunks = []
    t = 0
    for isi in isis:
        t += isi
        chunks.append(isi)
        if t > Tmax:
            ISIs_chunks.append(chunks.copy())
            chunks.clear()
            t = 0
    spectrum = []
    for f in fs:
        fws_real = []
        fws_img = []
        for isis in ISIs_chunks:
            fw = fourier_transformation_isis(f, isis)
            fws_real.append(fw.real)
            fws_img.append(fw.imag)
        spectrum.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))
    return spectrum


def inverse_gaussian(T, CV):
    ps = []
    ts = np.linspace(0, 2*T, 500)
    for t in ts:
        p = np.sqrt(T / (2 * np.pi * (CV**2) * (t ** 3))) * np.exp(-(t - T) ** 2 / (2 * T * (CV**2) * t))
        ps.append(p)
    return ts, ps