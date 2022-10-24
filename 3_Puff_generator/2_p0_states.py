import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax1])
    home = os.path.expanduser("~")

    ca_rest = 0.20
    ca = 0.20
    ip3 = 1.
    n = 5
    m = 3

    r0_opn = 0.1
    r_ref = 20.0
    r_cls = 50.
    r_opn_max = r0_opn*((1. + np.power(0.2, 3))/np.power(0.2, 3))*(2/1)
    r_opn = n * r_opn_max * (np.power(ca, 3) / (1. + np.power(ca, 3))) * (np.power(ip3, 3) / (1. + np.power(ip3, 3)))

    folder = home + "/Data/calcium_spikes_markov/ca_fix/"
    data = np.loadtxt(folder + f"puff_markov_cafix{ca:.2f}_ip{ip3:.2f}_tau1.00e+00_j1.00e+00_K1_5.dat")
    data_tmp = []
    for x in data:
        if x[2] == 0:
            data_tmp.append(x)
    data = data_tmp

    states = [i for i in range(n+m)]
    times_per_state = [[] for _ in range(n+m)]
    prob_per_state = [[] for _ in range(n+m)]

    for set1, set2 in zip(data[:-1], data[1:]): #set = [time, state, idx]
        state = int(set1[1])
        if state <= 0:
            #if the state is closed
            idx_state = state + (m-1)
        elif state > 0:
            #if state is open
            idx_state = -state
        time = set2[0] - set1[0] # calculate time spend in state x = set1[1]
        times_per_state[idx_state].append(time)

    time_per_state = [sum(times) for times in times_per_state]
    prob_per_state = [time/sum(time_per_state) for time in time_per_state]

    prob_per_state_theory = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)

    prob_per_state_theory2 = []
    for i in range(n+m):
        if i < m-1:
            p = 1/r_ref
        elif i == m-1:
            p = 1/r_opn
        else:
            p = (i-(m-1))/(n*r_cls)
        prob_per_state_theory2.append(p)
    norm = sum(prob_per_state_theory2)
    prob_per_state_theory2 = [p/norm for p in prob_per_state_theory2]

    ax1.plot(states, prob_per_state_theory2, lw=1, c=st.colors[2], zorder=1, label="Theory")
    ax1.scatter(states, prob_per_state, fc="w", ec=st.colors[0], s=20, zorder=2, label="Sim.")
    ax1.set_xticks(states)
    labels = [[] for i in range(n+m)]
    for i in range(n+m):
        if i < m:
            labels[i] = f"$0_{m-i}$"
        elif i >= m:
            labels[i] = n+m-i
    ax1.set_xticklabels(labels)
    ax1.set_yscale("log")
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$P_0(x)$")
    ax1.legend(fancybox=False, framealpha=1.0)
    plt.show()