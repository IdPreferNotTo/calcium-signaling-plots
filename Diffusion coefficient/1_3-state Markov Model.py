import numpy as np
from numpy.random import default_rng
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os

def steady_states(r1, r2, r3):
    denom = r1*r2 + r1*r3 + r2*r3
    p1 = r2*r3/denom
    p2 = r1*r3/denom
    p3 = r1*r2/denom
    return [p1, p2, p3]

def markov_model(ts, dt, r1, r2, r3):
    states = [] # time state
    state = 1
    for _ in ts:
        rng = np.random.uniform(0, 1)
        if state == 1:
            pt = r1
        elif state == 2:
            pt = r2
        elif state == 3:
            pt = r3
        if(rng < pt*dt):
            state += 1
            if state == 4:
                state = 1
        states.append(state)
    return states

def autocorrelation(xs, tstep, tmax):
    mean = np.mean(xs)
    xs = [x-mean for x in xs]
    correlations = []
    kmax = int(tmax/tstep)
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

def conditional_propability_over_t(x0, xtau, data):
    n=0
    for x in data:
        if x == xtau:
            n+=1
    p0_x0 = n/len(data)


    # get all indicies at which x0 is found
    idx_x0 = []
    for idx, x in enumerate(data):
        if x == x0:
          idx_x0.append(idx)

    p_xtau_after_x0_over_tau = []
    for k in range(200):
        # get all x values that are found a certain time tau=dt*(idx+k) after x0 was found
        x_after_x0 = []
        for idx in idx_x0:
            if idx+k < len(data)-1:
                x_after_x0.append(data[idx+k])

        # find the probability that xtau is found at time tau if x0 was found at time 0
        n_xtau_after_x0 = 0
        for x in x_after_x0:
            if x==xtau:
                n_xtau_after_x0 += 1
        p_xtau_after_x0 = n_xtau_after_x0/len(x_after_x0)
        p_xtau_after_x0_over_tau.append(p_xtau_after_x0)
    return p_xtau_after_x0_over_tau

if __name__ == "__main__":
    r1 = 0.1
    r2 = 0.2
    r3 = 0.3

    dt = 0.1
    ts = np.arange(0, 10000, step=dt)
    states = markov_model(ts, dt, r1, r2, r3)

    ps = [p1, p2, p3] = steady_states(r1, r2, r3)

    mean = 1*p1 + 2*p2 + 3*p3
    variance = (1-mean)**2*p1 + (2-mean)**2*p2 + (3-mean)**2*p3

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])

    tmax = 10
    taus, corr = autocorrelation(states, tstep=dt, tmax=tmax)
    ax.plot(taus, corr)
    ax.axhline(0, c="k")
    ax.set_xlim([0, tmax])

    condotional_prob = conditional_propability_over_t(1, 1, states)
    axins = inset_axes(ax, width="60%", height="50%", loc=1)
    axins.plot([dt*x for x in range(200)], condotional_prob)
    plt.show()