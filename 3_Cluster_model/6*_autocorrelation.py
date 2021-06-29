import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import gridspec

def delta(a, b):
    if a==b:
        return 1
    else:
        return 0


def steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n):
    A = np.array([[-n*r_ref, 0, 0, 0, 0, 0, 0, r_cls],
                  [n*r_ref, -n*r_ref, 0, 0, 0, 0, 0, 0],
                  [0, n*r_ref, -n*r_ref, 0, 0, 0, 0, 0],
                  [0, 0, n*r_ref, -n*r_opn, 0, 0, 0, 0],
                  [0, 0, 0, r_opn, -r_cls, 0, 0, 0],
                  [0, 0, 0, r_opn, r_cls, -r_cls, 0, 0],
                  [0, 0, 0, r_opn, 0, r_cls, -r_cls, 0],
                  [1, 1, 1, 1, 1, 1, 1, 1]])
    Ainv = np.linalg.inv(A)
    inhomgeneity = np.array([0, 0, 0, 0, 0, 0, 0, 1])
    p0s = Ainv.dot(inhomgeneity)
    return p0s


def f_from_k_invert_A(k, r_ref, r_opn, r_cls, m, n):
    A = np.array([[-n * r_ref, 0, 0, 0, 0, 0, 0, r_cls],
                  [n * r_ref, -n * r_ref, 0, 0, 0, 0, 0, 0],
                  [0, n * r_ref, -n * r_ref, 0, 0, 0, 0, 0],
                  [0, 0, n * r_ref, -n * r_opn, 0, 0, 0, 0],
                  [0, 0, 0, r_opn, -r_cls, 0, 0, 0],
                  [0, 0, 0, r_opn, r_cls, -r_cls, 0, 0],
                  [0, 0, 0, r_opn, 0, r_cls, -r_cls, 0],
                  [1, 1, 1, 1, 1, 1, 1, 1]])
    Ainv = np.linalg.inv(A)
    p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n)
    p0s[-1] = 0
    p0s = np.asarray(p0s)
    deltas = np.array([delta(k, 0), delta(k, 1), delta(k, 2), delta(k, 3), delta(k,4), delta(k, 5), delta(k, 6), 0])
    inhomgeneity = np.subtract(p0s, deltas)
    f_from_k = Ainv.dot(inhomgeneity)
    return f_from_k

def autocorrelation(data, dt):
    mean = np.mean(data)
    data = [x - mean for x in data]
    correlation = []
    ks = np.arange(0, 50)
    for k in ks:
        print(k)
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

def plot_correlation_function():
    ca_rest = 0.33
    home = os.path.expanduser("~")

    ca_fix1 = 0.3
    file = home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/puff_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca_fix1)
    data = np.loadtxt(file)
    t0 = data[0][0]
    data = [[set[0]-t0, set[1], set[2]] for set in data]

    # Preprocess data to contain data only from a single cluster
    data_single_cluster = []
    for set in data:
        if set[2] == 1:
            data_single_cluster.append(set)

    # Preprocess data so that refractory states are equal to 0
    data_single_cluster_without_ref = []
    for set in data_single_cluster:
        if set[1] > 0:
            data_single_cluster_without_ref.append(set)
        else:
            set[1] = 0
            data_single_cluster_without_ref.append(set)

    # Fill data
    data_single_cluster_filled = []
    n = 0
    dt = 0.01
    print("Fill Data")
    for t in np.linspace(0, 2_000, 200_000):
        if t < data_single_cluster_without_ref[n+1][0]:
            x = data_single_cluster_without_ref[n][1]
            data_single_cluster_filled.append(x)
        else:
            n += 1
            x = data_single_cluster_without_ref[n][1]
            data_single_cluster_filled.append(x)

    print("done")
    # Calculate correlaiton function
    taus, correlations1 = autocorrelation(data_single_cluster_filled, dt)
    integral_corr = 0
    for tau, corr in zip(taus, correlations1):
        integral_corr += corr*dt
    print(integral_corr)

    ca_fix2 = 0.9
    file = home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/puff_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca_fix2)
    data = np.loadtxt(file)
    t0 = data[0][0]
    data = [[set[0]-t0, set[1], set[2]] for set in data]

    # Preprocess data to contain data only from a single cluster
    data_single_cluster = []
    for set in data:
        if set[2] == 1:
            data_single_cluster.append(set)

    # Preprocess data so that refractory states are equal to 0
    data_single_cluster_without_ref = []
    for set in data_single_cluster:
        if set[1] > 0:
            data_single_cluster_without_ref.append(set)
        else:
            set[1] = 0
            data_single_cluster_without_ref.append(set)

    # Fill data
    data_single_cluster_filled = []
    n = 0
    dt = 0.01
    print("Fill Data")
    for t in np.linspace(0, 2_000, 200_000):
        if t < data_single_cluster_without_ref[n+1][0]:
            x = data_single_cluster_without_ref[n][1]
            data_single_cluster_filled.append(x)
        else:
            n += 1
            x = data_single_cluster_without_ref[n][1]
            data_single_cluster_filled.append(x)

    print("done")
    # Calculate correlaiton function
    taus, correlations2 = autocorrelation(data_single_cluster_filled, dt)


    r_opn = 0.13 * np.power(ca_fix1 / ca_rest, 3) * ((1 + ca_rest ** 3) / (1 + ca_fix1 ** 3))
    r_ref = 1.3 * np.power(ca_fix1 / ca_rest, 3) * ((1 + ca_rest ** 3) / (1 + ca_fix1 ** 3))
    r_cls = 50

    p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, 4, 4)

    xs = [0, 0, 0, 0, 4, 3, 2, 1]
    idxs = [0, 1, 2, 3, 4, 5, 6, 7]

    D_theory = 0
    for k in idxs:
        sum_over_i = 0
        f_from_k_to = f_from_k_invert_A(k, r_ref, r_opn, r_cls, 4, 4)
        for i in idxs:
            sum_over_i += xs[i] * f_from_k_to[i]
        D_theory += xs[k] * p0s[k] * sum_over_i
    print(D_theory)

    file_ca = home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/ca_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca_fix1)
    data = np.loadtxt(file_ca)
    ts, cas, jpuffs, adaps = np.transpose(data)
    var = 0.1 * np.var(jpuffs) /(2.*10)
    print(var)

    coarse_jpuff = []
    f = 10
    for i in range(0, len(jpuffs), f):
        average = sum(jpuffs[i:i + f]) / f
        coarse_jpuff.append(average)
    var2 = 1.0 * np.var(coarse_jpuff)/ (2. * 10)
    print(var2)

    for i in range(0, len(jpuffs), 100):
        average = sum(jpuffs[i:i + f]) / f
        coarse_jpuff.append(average)
    var3 = 10.0 * np.var(coarse_jpuff)/ (2. * 10)
    print(var3)

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 4))
    gs = gridspec.GridSpec(1, 2)

    ax1 = fig.add_subplot(gs[0])
    ax1.set_ylabel(r"$C_{xx}(\tau)$")
    ax1.set_xlabel(r"$\tau$")
    ax1.set_xlim([0, 0.5])
    ax1.axhline(0, ls=":", c="C7")
    ax1.plot(taus, correlations1, label="c = {:.1f}".format(ca_fix1))
    ax1.fill_between(taus, correlations1, 0, facecolor="C0", alpha=0.5, zorder=2)
    ax1.legend()

    ax2 = fig.add_subplot(gs[1])
    ax2.set_xlabel(r"$\tau$")
    ax2.set_xlim([0, 0.5])
    ax2.axhline(0, ls=":", c="C7")
    ax2.plot(taus, correlations2, label="c = {:.1f}".format(ca_fix2))
    ax2.fill_between(taus, correlations2, 0, facecolor="C0", alpha=0.5, zorder=2)
    ax2.legend()

    plt.savefig(home + "/Data/Calcium/Plots/cluster_correlation_ca{:.2f}_and_ca{:.2f}.pdf".format(ca_fix1, ca_fix2), transparent=True)
    plt.show()

if __name__ == "__main__":
    plot_correlation_function()