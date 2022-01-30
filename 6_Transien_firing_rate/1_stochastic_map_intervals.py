import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os
from scipy.optimize import curve_fit

import functions as fc
import styles as st

def quadratic_func(x, a, b, c):
    return (1/2)*a*(x**2) + b*x + c

def lin_func(x, a, b):
    return a*x + b

def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))


if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    taua = 829
    ampa = 0.0309
    tau = 10.5
    j = 0.0146
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spikes = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    data = np.loadtxt(folder + file)
    ISIs  = np.loadtxt(folder + file_spikes)
    ts, cas, jpuffs, adaps = np.transpose(data)
    spike_times = []
    cR = 0.33
    cT = 1.0
    Ai_sim = []

    for t, ca, jpuff, adap, adap_after in zip(ts[:-1], cas[:-1], jpuffs[:-1], adaps[:-1], adaps[1:]):
        if ca == 1.00:
            spike_times.append(t)
            Ai_sim.append(1. - adap_after)

    T_infty = np.mean(ISIs[100:])
    T0 = ISIs[0]
    a_infty = ampa/(1 - (1-ampa)*np.exp(-T_infty/taua))
    c_er_infty = 1 - a_infty
    # SOLVE 2D STOCHASTIC MAP

    ISImax = 50

    ais = []
    for i in range(100):
        x =  (1-ampa)*np.exp(-T_infty/taua)
        ai =  ampa*(1. - np.power(x, i))/(1. - x)
        ais.append(ai)

    taueff = taua/(T_infty - taua*np.log(1-ampa))
    Ti_theo = []
    for i in range(ISImax):
        Ti_theo.append(T_infty - (T_infty - T0)*np.exp(-i/taueff))

    Ais_theo = []
    for i in range(ISImax):
        Ais_theo.append(a_infty - (a_infty - 0)*np.exp(-i/taueff))

    ax1.plot(ts, adaps, lw=1, c=st.colors[1])
    ax1.axhline(c_er_infty, ls=":", lw=1, c="k")
    ax1.set_ylabel(r"$c_{ER}(t)$")
    ax1.set_xlabel("$t$")
    ax1.set_xlim([0, 1000])
    ax1.set_ylim([0.5, 1])

    ax3.plot(range(ISImax), Ais_theo, c="C3")
    ax3.scatter(range(ISImax), Ai_sim[:ISImax], s=10, fc="w", ec=st.colors[1], label="sim.")
    ax3.scatter(range(ISImax), ais[:ISImax], marker="^", s=10, fc="w", ec="C1", zorder=1, label="approx.")
    ax3.axhline(a_infty, ls=":", lw=1, c="k")
    ax3.text(10, a_infty, "$a_\infty$", va="bottom")
    ax3.set_xlabel("$i$")
    ax3.set_ylim(0, 0.5)
    ax3.set_ylabel("$a_i = 1 - c_{ER}(t_i)$")

    ax2.scatter(range(ISImax), ISIs[:ISImax], s=10, fc="w", ec=st.colors[1])
    #ax5.scatter(range(ISImax), Tis_theo[1:ISImax + 1], marker="^", s=10, fc="w", ec="C1", zorder=1)
    ax2.axhline(T_infty, lw = 1, ls=":", c="C7")
    ax2.axhline(T0, lw = 1, ls=":", c="C7")
    ax2.text(3, T_infty - 5, "$T_\infty$", va="top")
    ax2.text(10, T0, "$T_0$", va="bottom")
    ax2.set_xlabel("$i$")
    ax2.set_ylabel("$T_i$")
    ax2.arrow(18, T0, 0, T_infty - T0, fc="k", length_includes_head=True, head_width=1,
              head_length=1.0, lw=0.75, clip_on=False, zorder=5)
    ax2.arrow(18, T_infty, 0, -(T_infty - T0), fc="k", length_includes_head=True, head_width=1,
              head_length=1.0, lw=0.75, clip_on=False, zorder=5)
    ax2.text(20, (T0 + T_infty) / 2, r"$\Delta T $", ha="left")
    ax2.plot(range(ISImax), Ti_theo, c="C3")

    ax4.set_ylim([0, 100])
    ax4.set_xlabel("$a_\infty - a_i$")
    ax4.set_ylabel("$T_i$")
    ax4.scatter([a_infty - a for a in Ai_sim[:ISImax]], ISIs[:ISImax], s=10, fc="w", ec=st.colors[1])

    nr_isis = 20
    start_interpolation_at = 3
    index_isis = np.arange(nr_isis)
    popt, pcov = curve_fit(quadratic_func, [a_infty - a for a in Ai_sim[start_interpolation_at:ISImax]], ISIs[start_interpolation_at:ISImax], p0=(-1, -100, 100))

    print(popt)
    print(np.mean([aj - ai for ai, aj in zip(Ai_sim[0:20], Ai_sim[1:21])]))
    xs = np.linspace(0, a_infty - Ai_sim[0], 100)
    ys = []
    for x in xs:
        ys.append((1/2)*popt[0]*(x**2) + popt[1]*x + popt[2])
    ax4.plot(xs, ys, c="C3")
    print((1/2)*popt[0]/popt[1])
    plt.savefig(home + "/Data/Calcium/Plots/13_map_ISIs.pdf", transparent=True)
    plt.show()

    print(f"$T_\infty/T_0={T_infty/T0}$")

