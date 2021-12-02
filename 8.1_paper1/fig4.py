import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, cm
import os

import styles as st
import functions as fc

def get_mean_opn_channel(data):
    # data = [time, state, idx]
    data_single_cluster = []
    for set in data:
        if set[2] == 0:
            data_single_cluster.append(set)
    data_single_cluster_no_ref = []
    for set in data_single_cluster:
        if set[1] < 0:
            data_single_cluster_no_ref.append([set[0], 0, set[2]])
        else:
            data_single_cluster_no_ref.append(set)
    A = 0
    for set1, set2 in zip(data_single_cluster_no_ref[:-1], data_single_cluster_no_ref[1:]):
        time = set2[0] - set1[0]
        state = set1[1]
        A += time*state
    return A/1000


def get_var_open_channel(data):
    # data = [time, state, idx]
    data_single_cluster = []
    for set in data:
        if set[2] == 0:
            data_single_cluster.append(set)
    data_single_cluster_no_ref = []
    for set in data_single_cluster:
        if set[1] < 0:
            data_single_cluster_no_ref.append([set[0], 0, set[2]])
        else:
            data_single_cluster_no_ref.append(set)
    A2 = 0
    for set1, set2 in zip(data_single_cluster_no_ref[:-1], data_single_cluster_no_ref[1:]):
        time = set2[0] - set1[0]
        state = set1[1]
        A2 += time*(state**2)
    return A2/1000


def r_opn_full(x, y, n, a, b, c, K0, rmax):
    K = K0*(y**3 / (1 + y**3))
    r_opn_full = n * rmax * (np.power(x, a) / (1. + np.power(x, a))) * (np.power(K, c) / (np.power(K, c) + np.power(x, c)))
    return r_opn_full


def r_opn(x, y, n, a, b, rmax):
    r_opn = n * rmax * (np.power(x, a)/(1. + np.power(x, a)))*(np.power(y, b)/(1. + np.power(y, b)))
    return r_opn


def r_ref(x, y, n, a, b, rmax):
    r_ref = n * rmax * (np.power(x, a) / (1. + np.power(x, a))) * (np.power(y, b) / (1. + np.power(y, b)))
    return r_ref


if __name__ == "__main__":
    home = os.path.expanduser("~")
    n_cha = 5
    m_ref = 4
    a = 3
    b = 3
    c = 3
    K = 10
    r0_opn = 0.13
    r0_ref = 1.30
    r0_opn_max = r0_opn*((1. + np.power(0.33, 3))/np.power(0.33, 3))*(2/1)
    r0_ref_max = r0_ref*((1. + np.power(0.33, 3))/np.power(0.33, 3))*(2/1)

    cas = np.linspace(0.2 , 0.99, 80)
    mean_x0s = []
    mean_x0s_sim = []
    var_x0s = []
    var_x0s_sim = []
    var_jpuffs = []
    mean_jpuffs = []
    D_jpuffs = []
    for ca in cas:
        print(f"{ca:.2f}")
        folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/ca_fix/"
        file = f"puff_markov_cafix{ca:.2f}_ip1.00_tau1.00e+00_j1.00e+00_N10_0.dat"
        data = np.loadtxt(folder + file)
        mean_x0_sim = get_mean_opn_channel(data)
        var_x0_sim =  get_var_open_channel(data) - mean_x0_sim**2
        mean_x0s_sim.append(mean_x0_sim)
        var_x0s_sim.append(var_x0_sim)
        Popen = np.power(1 + 2 * np.power(n_cha * (n_cha + 1), -1, dtype=float) * np.power(0.33 / ca, 3) * (1 + ca ** 3) / (1 + 0.33 ** 3) * 10 / 0.02, -1)
        mean_x0 = ((n_cha + 2) / 3) * Popen
        mean_x0s.append(mean_x0)
        var_x0 = Popen*(n_cha + 1)*(n_cha + 2)/ 6 - mean_x0**2
        var_x0s.append(var_x0)

        file_ca = f"ca_markov_cafix{ca:.2f}_ip1.00_tau1.00e+00_j1.00e+00_N10_0.dat"
        data_ca = np.loadtxt(folder + file_ca)
        ts, cas_dyn, jpuffs, adap = np.transpose(data_ca)
        cgrain_factor2 = 50
        jpuffs_f2 = fc.coarse_grain_list(jpuffs, cgrain_factor2)
        mean_jpuff = np.mean(jpuffs_f2)/10
        mean_jpuffs.append(mean_jpuff)
        var_jpuff = cgrain_factor2*0.1*np.var(jpuffs_f2)/10
        var_jpuffs.append(var_jpuff)
        D1 = fc.intensity_puff_single(ca, 5, 4, 1)
        D_jpuffs.append(2*D1)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3*2.25, 2))
    grids = gridspec.GridSpec(1, 3)
    ax1 = fig.add_subplot(grids[0])
    ax2 = fig.add_subplot(grids[1])
    ax3 = fig.add_subplot(grids[2])
    axis = [ax1, ax2, ax3]
    st.remove_top_right_axis(axis)

    ax1.text(0.1, 0.95, "A", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, "B", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.1, 0.95, "C", fontsize=13, transform=ax3.transAxes, va='top')

    ax1.set_xlim([0.2, 1])
    ax1.set_xticks([0.33, 1])
    ax1.set_xticklabels(["$c_R$", "$c_T$"])
    ax1.set_xlabel("$c_i$")
    ax1.set_ylabel(r"$\langle x \rangle$")
    ax1.plot(cas, mean_x0s_sim, lw=1, color=st.colors[0])
    ax1.plot(cas, mean_x0s, lw=1, color=st.colors[2])

    ax2.set_xlim([0.2, 1])
    ax2.set_xticks([0.33, 1])
    ax2.set_xticklabels(["$c_R$", "$c_T$"])
    ax2.set_xlabel("$c_i$")
    ax2.set_ylabel(r"$\langle \Delta x^2 \rangle$")
    ax2.plot(cas, var_x0s_sim, lw=1, color=st.colors[0])
    ax2.plot(cas, var_x0s, lw=1, color=st.colors[2])

    ax3.set_xlim([0.2, 1])
    ax3.set_xticks([0.33, 1])
    ax3.set_xticklabels(["$c_R$", "$c_T$"])
    ax3.set_xlabel("$c_i$")
    ax3.set_ylabel(r"$J_{\rm puff}$")
    ax3.plot(cas, mean_jpuffs, lw=1, color=st.colors[0])
    ax3.plot(cas, mean_x0s, lw=1, color=st.colors[2])
    D_upper = [mean + np.sqrt(var) for mean, var in zip(mean_x0s, D_jpuffs)]
    D_lower = [mean - np.sqrt(var) for mean, var in zip(mean_x0s, D_jpuffs)]
    ax3.plot(cas, D_upper, lw=1, color=st.colors[2])
    ax3.plot(cas, D_lower, lw=1, color=st.colors[2])
    var_upper = [mean + np.sqrt(var) for mean, var in zip(mean_jpuffs, var_jpuffs)]
    var_lower = [mean - np.sqrt(var) for mean, var in zip(mean_jpuffs, var_jpuffs)]
    ax3.fill_between(cas, var_upper, var_lower, alpha=0.5, color=st.colors[0])
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig4.pdf",transparent=True)
    plt.show()