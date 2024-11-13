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
    h = 1.5 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    st.remove_top_right_axis([ax1, ax2])
    ax1.text(0.1, 0.95, r"A", fontsize=11, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, r"B", fontsize=11, transform=ax2.transAxes, va='top')

    #ax1a = ax1.twinx()
    #ax2a = ax2.twinx()

    K = 5
    N = 5
    tau = 1
    mu = fc.mean_jp_single_theory(0.5, N, 3, 1)
    tauer = 300
    eps = 0.03
    ds = [0.05] #[-0.1, -0.05, 0.0, 0.05, 0.1]
    colors = []
    bts = np.linspace(0, 20, 101)
    for i, d in enumerate(ds):
        p = (1. + d) * (0.5 - 0.2) / (K * mu * tau)
        Ts_buffer = []
        CVs_buffer = []
        for bt in bts:
            isis_buffer = np.loadtxt(
            home + f"/Data/calcium/markov/adaptive/fast_buffer/"
                   f"spike_times_markov_buffer_bt{bt:.2f}_ip1.00_taua{tauer:.2e}_ampa{eps:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat")
            mean = np.mean(isis_buffer)
            std = np.std(isis_buffer)
            Ts_buffer.append(mean)
            CVs_buffer.append(std/mean)

        ax2.plot(bts, CVs_buffer, c=st.colors[1], zorder=2)
        ax1.plot(bts, Ts_buffer, c=st.colors[1], zorder=2, label="mean-driven")

    p = 0.95 * (0.5 - 0.2) / (K * mu * tau)
    print(p)
    Ts_buffer = []
    CVs_buffer = []
    for bt in bts:
        isis_buffer = np.loadtxt(
        home + f"/Data/calcium/markov/adaptive/fast_buffer/"
               f"spike_times_markov_buffer_bt{bt:.2f}_ip1.00_taua{tauer:.2e}_ampa{eps:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat")
        mean = np.mean(isis_buffer)
        std = np.std(isis_buffer)
        Ts_buffer.append(mean)
        CVs_buffer.append(std/mean)
    ax2.plot(bts, CVs_buffer, c="#009957ff", zorder=2)
    ax1.plot(bts, Ts_buffer, c="#009957ff", zorder=2, label="excitable")

    #ax1.set_title("mean-driven")
    ax1.set_ylim([0.0, 1_000])
    ax1.set_xlabel(r"$b_T$")
    ax1.set_ylabel(r"$\langle T \rangle$ / s")

    #ax2.set_title("excitable")
    ax2.set_ylim([0.2, 0.7])
    ax2.set_xlabel(r"$b_T$")
    ax2.set_ylabel(r"$C_V$")

    #ax1.axvspan(3.2, 13, alpha=0.2, color="C7")
    #ax2.axvspan(3.2, 13, alpha=0.2, color="C7")

    leg = ax1.legend(fancybox=False, fontsize=8, edgecolor="k", bbox_to_anchor=(0.0, 1.1, 2.45, .0), loc=3,
                      ncol=2, mode="expand", borderaxespad=0)
    leg.get_frame().set_linewidth(.75)
    #ax1a.set_ylim([0, 300])
    #ax2a.set_ylim([0, 1000])

    #secax1 = ax1.secondary_yaxis("top", functions=(tau_to_bt, bt_to_tau))
    #secax1.set_xlabel("$b_T$")

    #secax2 = ax2.secondary_xaxis("top", functions=(tau_to_bt, bt_to_tau))
    #secax2.set_xlabel("$b_T$")

    home = os.path.expanduser("~")
    plt.savefig(home +"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/SUB2/figures/fig13.pdf")
    plt.show()