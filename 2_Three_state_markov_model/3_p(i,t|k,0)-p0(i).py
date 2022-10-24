import numpy as np
from numpy.random import default_rng
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os


def steady_states_theory(r1, r2, r3):
    denom = r1 * r2 + r1 * r3 + r2 * r3
    p1 = r2 * r3 / denom
    p2 = r1 * r3 / denom
    p3 = r1 * r2 / denom
    return [p1, p2, p3]

def steady_state_sim(data):
    N = 0
    n1 = 0
    n2 = 0
    n3 = 0
    for x in data:
        N += 1
        if x == 2:
            n1 += 1
        if x == 3:
            n2 += 1
        if x == 1:
            n3 += 1
    p1 = n1/N
    p2 = n2/N
    p3 = n3/N
    return [p1, p2, p3]

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

def conditional_propability_over_t(k, i, tau, dt, data):
    xs = [2, 3, 1]
    # calculate the steady state probability of a xk
    n = 0
    N = 0
    for x in data:
        N += 1
        if x == xs[i-1]:
            n += 1
    p0_xi = n / N

    # get all indices at which x0 is found
    idx_k = []
    for idx, x in enumerate(data):
        if x == xs[k-1]:
            idx_k.append(idx)

    smax = int(tau / dt)
    prob_i_at_t_given_k_at_0 = []
    for s in range(smax):
        # get all x values that are found a certain time tau=dt*(idx+k) after xk was found
        xs_after_xk = []
        for idx in idx_k:
            if idx + s < len(data) - 1:
                xs_after_xk.append(data[idx + s])

        # find the probability that xtau is found at time tau if x0 was found at time 0
        n_i_after_k = 0
        for x in xs_after_xk:
            if x == xs[i-1]:
                n_i_after_k += 1
        p_i_after_k = n_i_after_k / len(xs_after_xk) - p0_xi
        prob_i_at_t_given_k_at_0.append(p_i_after_k)
    return prob_i_at_t_given_k_at_0

def F_particular_from_k_to_i(k, i, r1, r2, r3):
    def delta(i, j):
        if i==j:
            return 1
        else:
            return 0
    [p1, p2, p3] = steady_states_theory(r1, r2, r3)
    if i == 1:
        return (-p1 + delta(1, k))/r1
    if i == 2:
        return (p3 - delta(3, k))/r2
    if i == 3:
        return 0

if __name__ == "__main__":
    home = os.path.expanduser("~")
    r1 = 0.1
    r2 = 0.2
    r3 = 0.3

    dt = 0.1
    print("Loading data...")
    data = np.loadtxt(home + "/Data/3_state_Markov_model/3_state_t1e+05.dat")
    ts, states = np.transpose(data)
    print("done")
    ps = [p1, p2, p3] = steady_states_theory(r1, r2, r3)
    # r = normalization of the steady state probability
    r = r1*r2*r3/(r1*r2 + r1*r3 + r2*r3)

    k = 1
    # Fpart is the particular solution of the inhomogeneous system of linear equations
    Fpart_from_k_to_1 = F_particular_from_k_to_i(k, 1, r1, r2, r3)
    Fpart_from_k_to_2 = F_particular_from_k_to_i(k, 2, r1, r2, r3)
    Fpart_from_k_to_3 = F_particular_from_k_to_i(k, 3, r1, r2, r3)
    sumFpart = (Fpart_from_k_to_1 + Fpart_from_k_to_2 + Fpart_from_k_to_3)
    print(sumFpart)

    # F21 is the integrate probability from 1 to 2.
    F_from_k_to_1 = Fpart_from_k_to_1 - p1 * sumFpart
    F_from_k_to_2 = Fpart_from_k_to_2 - p2 * sumFpart
    F_from_k_to_3 = Fpart_from_k_to_3 - p3 * sumFpart
    print(F_from_k_to_1, F_from_k_to_2, F_from_k_to_3)

    conditional_prob_from_k_to_1 = conditional_propability_over_t(k, 1, tau=20, dt=dt, data=states)
    print("1 to 1 - done")
    print(sum(conditional_prob_from_k_to_1)*dt)
    conditional_prob_from_k_to_2 = conditional_propability_over_t(k, 2, tau=20, dt=dt, data=states)
    print("1 to 2 - done")
    conditional_prob_from_k_to_3 = conditional_propability_over_t(k, 3, tau=20, dt=dt, data=states)
    print("1 to 3 - done")
    conditional_prob_k_to_any = [x+y+z for (x, y, z) in zip(conditional_prob_from_k_to_1, conditional_prob_from_k_to_2, conditional_prob_from_k_to_3)]
    integral_prob_1_to_1 = sum(conditional_prob_from_k_to_1)*dt
    integral_prob_1_to_2 = sum(conditional_prob_from_k_to_2)*dt
    integral_prob_1_to_3 = sum(conditional_prob_from_k_to_3)*dt

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.plot([dt * x for x in range(200)], conditional_prob_from_k_to_1, label="$k = 1$, $i = 1$, $F(1, 1) = {:.2f}/{:.2f}$".format(integral_prob_1_to_1, F_from_k_to_1))
    ax.plot([dt * x for x in range(200)], conditional_prob_from_k_to_2, label="$k = 1$, $i = 2$, $F(1, 2) = {:.2f}/{:.2f}$".format(integral_prob_1_to_2, F_from_k_to_2))
    ax.plot([dt * x for x in range(200)], conditional_prob_from_k_to_3, label="$k = 1$, $i = 3$, $F(1, 3) = {:.2f}/{:.2f}$".format(integral_prob_1_to_3, F_from_k_to_3))
    ax.plot([dt * x for x in range(200)], conditional_prob_k_to_any, c="k")
    ax.set_xlim([0, 20])
    ax.legend()
    ax.set_ylabel(r"$p(x_i, \tau | x_k, 0) - p(x_i)$")
    ax.set_xlabel(r"$\tau$")
    #plt.savefig(home + "/Data/Calcium/Plots/puff_model_transition_probability_{:d}_to_k.pdf".format(x0), transparent=True)
    plt.show()
