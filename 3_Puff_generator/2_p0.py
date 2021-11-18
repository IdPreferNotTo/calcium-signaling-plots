import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax1])
    home = os.path.expanduser("~")

    ca_rest = 0.33
    ca = 0.4
    ip3 = 1.
    n = 5
    m = 4

    meanI = 10 # for a single channelCaFix
    r_opn = 0.13 * np.power(ca / ca_rest, 3) * (1 + ca_rest ** 3) / (1 + ca ** 3)
    r_ref = 1.3 * np.power(ca / ca_rest, 3) * (1 + ca_rest ** 3) / (1 + ca ** 3)
    r_cls = 50
    pnorm = ((m-1) * r_cls * r_opn + r_cls * r_ref + (n + 1) * (n / 2) * r_opn * r_ref) / (r_cls * r_opn * r_ref)

    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/"
    data = np.loadtxt(folder + f"puff_markov_cafix{ca:.2f}_ip{ip3:.2f}_tau1.00e+00_j1.00e+00_N10_0.dat")
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

    prob_per_state_theory = []
    for i in range(1, n+m+1):
        # i in [1, n+m] including n+m
        if i < m:
            # i < m refractory states 0_m - 0_2
            prob_per_state_theory.append(1/r_ref)
        elif i == m:
            # i == m closed state 0_1
            prob_per_state_theory.append(1/r_opn)
        elif i > m:
            # i > m open states with n open channel
            n = i - m
            prob_per_state_theory.append(n/r_cls)
    # prob still has to be normalized
    prob_total = sum(prob_per_state_theory)
    prob_per_state_theory = [p/prob_total for p in prob_per_state_theory]
    p_open = sum(prob_per_state_theory)

    ax1.plot(states, prob_per_state_theory, lw=1, c=st.colors[2], zorder=1, label="Theory")
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
    plt.savefig(home + "/Data/Calcium/Plots/puff_gen_steady_state_prob.pdf".format(ca))
    plt.show()