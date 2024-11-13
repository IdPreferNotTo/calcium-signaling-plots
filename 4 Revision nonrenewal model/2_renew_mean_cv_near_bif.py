import numpy as np
import matplotlib.pyplot as plt
import os

import functions as fc
import default_parameters as df
import styles as st

mu = fc.mean_jp_single_theory(0.5, 5, 3, 1)
print(mu)

def tau_to_bt(xs):
    return 5.*(xs/1. -1.) #K* (tau/tau* - 1.)

def bt_to_tau(xs):
    return (1+ xs/5.)*1. # (1 + bT/K*)tau*


if __name__ == "__main__":
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    width, height = st.figsize([3, 2], "large")
    fig, axes = plt.subplots(nrows=2, ncols=3, layout="constrained", figsize=(width, height))
    ax1 = axes[0, 2]
    ax2 = axes[1, 2]
    ax3 = axes[0, 1]
    ax4 = axes[1, 1]
    ax5 = axes[0, 0]
    ax6 = axes[1, 0]
    ax1.set_title("Below bif. - Excitable")
    ax3.set_title("On bif.")
    ax5.set_title("Above bif. - Mean-driven")
    K = 5
    N = 5
    taus = np.linspace(1, 10, num=10)
    mu = fc.mean_jp_single_theory(0.5, N, 3, 1)

    ps = (0.5 - 0.2) / (K * mu * 1.05*taus)
    T0s = []
    CVs = []
    for tau, p in zip(taus, ps):
        isis = df.load_spike_times_markov(tau, p, K=K, N=N, cer=False)
        mean = np.mean(isis)
        std = np.std(isis)
        T0s.append(mean)
        CVs.append(std/mean)
    ax1.plot(taus, T0s, c=st.colors[0])
    ax2.plot(taus, CVs, c=st.colors[0])

    T0s_buffer = []
    CVs_buffer = []
    bts = [0, 10,20, 30, 40]
    tau_ref = taus[0]
    p_ref = ps[0]
    for bt in bts:
        isis_buffer = np.loadtxt(
        home + f"/Data/calcium/markov/renewal/fast_buffer/"
               f"spike_times_markov_buffer_bt{bt:.2f}_ip1.00_tau{tau_ref:.2e}_j{p_ref:.2e}_K{K:d}_{N:d}.dat")
        mean = np.mean(isis_buffer)
        std = np.std(isis_buffer)
        T0s_buffer.append(mean)
        CVs_buffer.append(std/mean)
    ax1.scatter([1. + bt/5. for bt in bts], T0s_buffer, fc="w", ec=st.colors[0], zorder=2)
    ax2.scatter([1. + bt/5. for bt in bts], CVs_buffer, fc="w", ec=st.colors[0], zorder=2)

    # ps = (0.5 - 0.2) / (K * mu * taus)
    # T0s = []
    # CVs = []
    # for tau, p in zip(taus, ps):
    #     isis = df.load_spike_times_markov(tau, p, K=K, N=N, cer=False)
    #     mean = np.mean(isis)
    #     std = np.std(isis)
    #     T0s.append(mean)
    #     CVs.append(std/mean)
    # ax3.plot(taus, T0s, c=st.colors[0])
    # ax4.plot(taus, CVs, c=st.colors[0])

    ps = (0.5 - 0.2) / (K * mu * 0.95*taus)
    print(ps)
    T0s = []
    CVs = []
    for tau, p in zip(taus, ps):
        isis = df.load_spike_times_markov(tau, p, K=K, N=N, cer=False)
        mean = np.mean(isis)
        std = np.std(isis)
        T0s.append(mean)
        CVs.append(std/mean)
    ax5.plot(taus, T0s, c=st.colors[0])
    ax6.plot(taus, CVs, c=st.colors[0])

    T0s_buffer = []
    CVs_buffer = []
    bts = [0, 10,20, 30, 40]
    p_ref = ps[0]
    for bt in bts:
        isis_buffer = np.loadtxt(
        home + f"/Data/calcium/markov/renewal/fast_buffer/"
               f"spike_times_markov_buffer_bt{bt:.2f}_ip1.00_tau{tau_ref:.2e}_j{p_ref:.2e}_K{K:d}_{N:d}.dat")
        mean = np.mean(isis_buffer)
        std = np.std(isis_buffer)
        T0s_buffer.append(mean)
        CVs_buffer.append(std/mean)
    ax5.scatter([1. + bt/5. for bt in bts], T0s_buffer, fc="w", ec=st.colors[0], zorder=2)
    ax6.scatter([1. + bt/5. for bt in bts], CVs_buffer, fc="w", ec=st.colors[0], zorder=2)

    ax1.set_xlabel(r"$\tau$")
    ax1.set_ylabel(r"$T_0$")
    ax2.set_xlabel(r"$\tau$")
    ax2.set_ylabel(r"$C_V(T_0)$")

    ax3.set_xlabel(r"$\tau$")
    ax3.set_ylabel(r"$T_0$")
    ax4.set_xlabel(r"$\tau$")
    ax4.set_ylabel(r"$C_V(T_0)$")

    ax5.set_xlabel(r"$\tau$")
    ax5.set_ylabel(r"$T_0$")
    ax6.set_xlabel(r"$\tau$")
    ax6.set_ylabel(r"$C_V(T_0)$")

    secax1 = ax1.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    secax2 = ax2.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    secax1.set_xlabel("$b_T$")
    secax2.set_xlabel("$b_T$")

    secax5 = ax5.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    secax6 = ax6.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    secax5.set_xlabel("$b_T$")
    secax6.set_xlabel("$b_T$")
    home = os.path.expanduser("~")
    plt.savefig(home + f"/Desktop/hek_fit/renew_mean_cv_close_to_bifurcation.png", dpi=300)
    plt.show()