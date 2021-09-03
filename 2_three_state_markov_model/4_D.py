import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os


def steady_states_theory(r1, r2, r3):
    denom = r1 * r2 + r1 * r3 + r2 * r3
    p1 = r2 * r3 / denom
    p2 = r1 * r3 / denom
    p3 = r1 * r2 / denom
    return [p1, p2, p3]

def f_particular_from_k_to_i(k, i, r1, r2, r3):
    def delta(a, b):
        if a==b:
            return 1
        else:
            return 0
    [p1, p2, p3] = steady_states_theory(r1, r2, r3)
    if i == 1:
        return (-p1 + delta(1, k))/r1
    if i == 2:
        return (-p1 + delta(1, k) - p2 + delta(2, k))/r2
    if i == 3:
        return 0

def autocorrelation(xs, tstep, tmax):
    mean = np.mean(xs)
    xs = [x - mean for x in xs]
    correlations = []
    kmax = int(tmax / tstep)
    ks = np.arange(0, kmax)
    for k in ks:
        print(k)
        corr = 0
        if k == 0:
            for x, y in zip(xs, xs):
                corr += x * y
        else:
            for x, y in zip(xs[:-k], xs[k:]):
                corr += x * y
        correlations.append(corr / len(xs[k:]))
    return [[tstep * k for k in ks], correlations]

def F_from_xk_to_xi(k, i, p0_xi):
    idxs = [1,2,3]
    sum_f_part = 0
    for l in idxs:
        sum_f_part += f_particular_from_k_to_i(k, l, r1, r2, r3)
    f_part_from_k_to_i = f_particular_from_k_to_i(k, i, r1, r2, r3)
    f_from_k_to_i = f_part_from_k_to_i -  p0_xi*sum_f_part
    return f_from_k_to_i

if __name__ == "__main__":
    home = os.path.expanduser("~")
    r1 = 0.1
    r2 = 0.2
    r3 = 0.3

    data = np.loadtxt(home + "/Data/3_state_Markov_model/3_state_t1e+05.dat")
    ts, states = np.transpose(data)

    ps = [p1, p2, p3] = steady_states_theory(r1, r2, r3)

    x0=0
    idxs = [1,2,3]
    xs = [2,3,1]

    mean = sum([x*p for x, p in zip(xs, ps)])
    D_theory = 0
    for k in idxs:
        sum_over_i = 0
        for i in idxs:
            sum_over_i += (xs[i-1] - mean)*f_particular_from_k_to_i(k, i, r1, r2, r3)
        D_theory += xs[k-1] * ps[k-1]*sum_over_i

    tmax = 20
    dt = 0.1
    taus, corr = autocorrelation(states, tstep=dt, tmax=tmax)

    D_sim = 0
    for c in corr:
        D_sim += c
    D_sim = D_sim*dt

    print(D_sim, D_theory)

    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.plot(taus, corr)
    ax.axhline(0, c="k")
    ax.set_xlim([0, tmax])
    ax.fill_between(taus, corr, 0, facecolor="C0", alpha=0.5, zorder=2, label="D = {:.2f} / {:.2f}".format(D_sim, D_theory))
    plt.legend()
    plt.savefig(home + "/Data/Calcium/Plots/Diffusion_coefficient_3_state_Markov_model.pdf".format(x0), transparent=True)
    plt.show()
