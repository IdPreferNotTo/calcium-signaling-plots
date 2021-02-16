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
        if x == 1:
            n1 += 1
        if x == 2:
            n2 += 1
        if x == 3:
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

def conditional_propability_over_t(x0, xtau, tau, dt, data):
    # calculate the steady state probability of a x0
    n = 0
    N = 0
    for x in data:
        N += 1
        if x == xtau:
            n += 1
    p0_xtau = n / N

    # get all indices at which x0 is found
    idx_x0 = []
    for idx, x in enumerate(data):
        if x == x0:
            idx_x0.append(idx)

    kmax = int(tau / dt)
    p_xtau_after_x0_over_tau = []
    for k in range(kmax):
        # get all x values that are found a certain time tau=dt*(idx+k) after x0 was found
        x_after_x0 = []
        for idx in idx_x0:
            if idx + k < len(data) - 1:
                x_after_x0.append(data[idx + k])

        # find the probability that xtau is found at time tau if x0 was found at time 0
        n_xtau_after_x0 = 0
        for x in x_after_x0:
            if x == xtau:
                n_xtau_after_x0 += 1
        p_xtau_after_x0 = n_xtau_after_x0 / len(x_after_x0) - p0_xtau
        p_xtau_after_x0_over_tau.append(p_xtau_after_x0)
    return p_xtau_after_x0_over_tau


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

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])

    # r = normalization of the steady state probability
    r = r1*r2*r3/(r1*r2 + r1*r3 + r2*r3)

    # Fpart is the particular solution of the inhomogeneous system of linear equations
    Fpart11 = (-p1 + 1)/r1
    Fpart21 = p3/r2
    Fpart31 = 0
    sumFpart = (Fpart11 + Fpart21 + Fpart31)

    # F21 is the integrate probability from 1 to 2.
    F11 = Fpart11 - p1 * sumFpart
    F21 = Fpart21 - p2 * sumFpart
    F31 = Fpart31 - p3 * sumFpart
    print(F11, F21, F31)

    conditional_prob_1_to_1 = conditional_propability_over_t(x0=1, xtau=1, tau=20, dt=dt, data=states)
    print("1 to 1 - done")
    conditional_prob_1_to_2 = conditional_propability_over_t(x0=1, xtau=2, tau=20, dt=dt, data=states)
    print("1 to 2 - done")
    conditional_prob_1_to_3 = conditional_propability_over_t(x0=1, xtau=3, tau=20, dt=dt, data=states)
    print("1 to 3 - done")
    conditional_prob_1_to_any = [x+y+z for (x, y, z) in zip(conditional_prob_1_to_1, conditional_prob_1_to_2, conditional_prob_1_to_3)]
    integral_prob_1_to_1 = sum(conditional_prob_1_to_1)*dt
    integral_prob_1_to_2 = sum(conditional_prob_1_to_2)*dt
    integral_prob_1_to_3 = sum(conditional_prob_1_to_3)*dt
    ax.plot([dt * x for x in range(200)], conditional_prob_1_to_1, label="$x_k = 1$, $x_i = 1$, $F(1, 1) = {:.2f}/{:.2f}$".format(integral_prob_1_to_1, F11))
    ax.plot([dt * x for x in range(200)], conditional_prob_1_to_2, label="$x_k = 1$, $x_i = 2$, $F(2, 1) = {:.2f}/{:.2f}$".format(integral_prob_1_to_2, F21))
    ax.plot([dt * x for x in range(200)], conditional_prob_1_to_3, label="$x_k = 1$, $x_i = 2$, $F(3, 1) = {:.2f}/{:.2f}$".format(integral_prob_1_to_3, F31))
    ax.plot([dt * x for x in range(200)], conditional_prob_1_to_any, c="k")
    ax.set_xlim([0, 20])
    ax.legend()
    ax.set_ylabel(r"$p(x_i, \tau | x_0, 0) - p(x_i)$")
    ax.set_xlabel(r"$\tau$")
    plt.savefig(home + "/Data/Calcium/Plots/puff_model_transition_probability_1_to_k.pdf", transparent=True)

    plt.show()
