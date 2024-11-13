import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import gridspec

import functions as fc
import default_parameters as df
import styles as st

mu = fc.mean_jp_single_theory(0.5, 5, 3, 1)
print(mu)

def tau_to_bt(xs):
    return 1.*(xs/1. -1.) #K* (tau/tau* - 1.)

def bt_to_tau(xs):
    return (1+ xs/1.)*1. # (1 + bT/K*)tau*


if __name__ == "__main__":
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.66 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1])

    ax1.text(0.05, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax3.text(0.05, 0.95, r"B", fontsize=11, transform=ax3.transAxes, va='top')

    #ax3.set_title("Excitable")
    #ax1.set_title("Mean-driven")
    K = 5
    N = 5
    taus = np.linspace(1, 10, num=50)
    mu = fc.mean_jp_single_theory(0.5, N, 3, 1)
    tauer = 300
    eps = 0.03
    ps = (0.5 - 0.2) / (K * mu * 1.05*taus)
    tau_ref = taus[0]
    p_ref = ps[0]
    T0s = []
    CVs = []
    for tau, p in zip(taus, ps):
        isis = df.load_spike_times_markov(tau, p, cer=True, taua=tauer, ampa=eps, K=K, N=N)
        mean = np.mean(isis)
        std = np.std(isis)
        T0s.append(mean)
        CVs.append(std/mean)
    ax3.plot(taus, T0s, c=st.colors[3])
    ax4.plot(taus, CVs, c=st.colors[3])

    T0s_buffer = []
    CVs_buffer = []
    bts = [0, 4, 8, 12, 16]
    for bt in bts:
        isis_buffer = np.loadtxt(
        home + f"/Data/calcium/markov/adaptive/fast_buffer/"
               f"spike_times_markov_buffer_bt{bt:.2f}_ip1.00_taua{tauer:.2e}_ampa{eps:.2e}_tau{tau_ref:.2e}_j{p_ref:.2e}_K{K:d}_{N:d}.dat")
        mean = np.mean(isis_buffer)
        std = np.std(isis_buffer)
        T0s_buffer.append(mean)
        CVs_buffer.append(std/mean)
    ax3.scatter([1. + bt / 2.25 for bt in bts], T0s_buffer, fc="w", ec=st.colors[3], zorder=2)
    ax4.scatter([1. + bt / 2.25 for bt in bts], CVs_buffer, fc="w", ec=st.colors[3], zorder=2)

    # ps = (0.5 - 0.2) / (K * mu * taus)
    # T0s = []
    # CVs = []
    # for tau, p in zip(taus, ps):
    #     isis = df.load_spike_times_markov(tau, p, cer=True, taua=tauer, ampa=eps, K=K, N=N)
    #     mean = np.mean(isis)
    #     std = np.std(isis)
    #     T0s.append(mean)
    #     CVs.append(std/mean)
    # ax3.plot(taus, T0s, c=st.colors[0])
    # ax4.plot(taus, CVs, c=st.colors[0])

    ps = (0.5 - 0.2) / (K * mu * 0.95 * taus)
    p_ref = ps[0]
    T0s = []
    CVs = []
    for tau, p in zip(taus, ps):
        isis = df.load_spike_times_markov(tau, p, cer=True, taua=tauer, ampa=eps, K=K, N=N)
        mean = np.mean(isis)
        std = np.std(isis)
        T0s.append(mean)
        CVs.append(std/mean)
    ax1.plot(taus, T0s, c=st.colors[1])
    ax2.plot(taus, CVs, c=st.colors[1])

    T0s_buffer = []
    CVs_buffer = []
    bts = [0, 4, 8, 12, 16]
    for bt in bts:
        isis_buffer = np.loadtxt(
        home + f"/Data/calcium/markov/adaptive/fast_buffer/"
               f"spike_times_markov_buffer_bt{bt:.2f}_ip1.00_taua{tauer:.2e}_ampa{eps:.2e}_tau{tau_ref:.2e}_j{p_ref:.2e}_K{K:d}_{N:d}.dat")
        mean = np.mean(isis_buffer)
        std = np.std(isis_buffer)
        T0s_buffer.append(mean)
        CVs_buffer.append(std/mean)
    ax1.scatter([1. + bt/2.25 for bt in bts], T0s_buffer, fc="w", ec=st.colors[1], zorder=2)
    ax2.scatter([1. + bt/2.25 for bt in bts], CVs_buffer, fc="w", ec=st.colors[1], zorder=2)

    #ax3.set_xlabel(r"$\tau$")
    ax3.set_xticks([0])
    #ax3.set_ylim([0, 1])
    ax4.set_ylim([0.2, 0.7])
    #ax4.set_xlim([0, 11])
    ax3.set_ylabel(r"$T$")
    ax4.set_xlabel(r"$\tau$")
    ax4.set_ylabel(r"$C_V(T)$")

    #ax1.set_xlabel(r"$\tau$")
    ax1.set_xticks([0])
    #ax1.set_ylim([])
    ax2.set_ylim([0.2, 0.7])
    #ax2.set_xlim([0, 11])
    ax1.set_ylabel(r"$T$")
    ax2.set_xlabel(r"$\tau$")
    ax2.set_ylabel(r"$C_V(T)$")

    secax1 = ax1.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    #secax2 = ax2.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    secax1.set_xlabel("$b_T$")
    #secax2.set_xlabel("$b_T$")

    secax3 = ax3.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    #secax4 = ax4.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    secax3.set_xlabel("$b_T$")
    #secax4.set_xlabel("$b_T$")


    home = os.path.expanduser("~")
    plt.savefig(home + f"/Desktop/hek_fit/refrac_mean_cv_close_to_bifurcation.pdf", dpi=300)
    plt.show()