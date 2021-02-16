import numpy as np
import time
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


if __name__ == "__main__":
    home = os.path.expanduser("~")
    r1 = 0.1
    r2 = 0.2
    r3 = 0.3

    data = np.loadtxt(home + "/Data/3_state_Markov_model/3_state.dat")
    ts, states = np.transpose(data)

    ps_theory = [p1_t, p2_t, p3_t] = steady_states_theory(r1, r2, r3)
    ps_sim = [p1_s, p2_s, p3_s] = steady_state_sim(states)

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.plot(ps_theory)
    ax.plot(ps_sim)
    plt.show()