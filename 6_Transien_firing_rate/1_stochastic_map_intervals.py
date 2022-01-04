import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os
from scipy.optimize import curve_fit, fsolve

import functions as fc
import styles as st
def lin_func(x, m):
    return m*x

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, (3/2)*6 / 2))
    gs = gridspec.GridSpec(3, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[0, 1])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[2, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4, ax5, ax6])

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
    t_isi = []
    ca_isi = []
    jpuff_isi = []
    adap_isi= []
    Ai = [1]
    for t, ca, jpuff, adap, adap_after in zip(ts[:-1], cas[:-1], jpuffs[:-1], adaps[:-1], adaps[1:]):
        t_isi.append(t)
        ca_isi.append(ca)
        jpuff_isi.append(jpuff)
        adap_isi.append(adap)
        if ca == 1.00:
            spike_times.append(t)
            Ai.append(1. - adap_after)
            ax1.plot(t_isi, adap_isi, lw=1, c="k")
            ax1.plot([t, t], [adap, adap_after], lw=1, c="C7")
            t_isi.clear()
            ca_isi.clear()
            jpuff_isi.clear()
            adap_isi.clear()
    T_infty = np.mean(ISIs[100:])
    T0 = ISIs[0]
    a_infty =  1/(1 + (ampa*taua/T_infty))
    b_infty = ampa/(1 - (1-ampa)*np.exp(-T_infty/taua))

    # SOLVE 2D STOCHASTIC MAP
    def func_Ti(Ti, T0, bi, taua):
        return Ti - T0 - bi*taua*(1-np.exp(-Ti/taua))

    bis_theo = [0]
    Tis_theo = [T0]
    for i in range(100):
        bi = bis_theo[-1]
        Ti = Tis_theo[-1]
        bj = ampa*(1 - np.power(np.exp(-T_infty/taua)*(1-ampa), i))/(1 - np.exp(-T_infty/taua)*(1-ampa))#bi*np.exp(-Ti/taua) + ampa*(1-bi*np.exp(-Ti/taua))
        Tj = T0 / (1 - bi)
        bis_theo.append(bj)
        Tis_theo.append(Tj)

    ISImax = 25
    ax1.axhline(a_infty, ls=":", lw=1, c="k")
    ax1.set_ylabel(r"$a(t)$")
    ax1.set_xlim([0, 1000])
    ax1.set_ylim([0.5, 1])

    ax4.scatter(range(ISImax), Ai[:ISImax], s=10, fc="w", ec="k")
    ax4.scatter(range(ISImax), bis_theo[:ISImax], marker="^", s=10, fc="w", ec="C1", zorder=1)
    ax4.axhline(b_infty, ls=":", lw=1, c="k")
    ax4.text(10, b_infty, "$b_\infty$", va="bottom")
    ax4.set_xlabel("$i$")
    ax4.set_ylim(0, 0.5)
    ax4.set_ylabel("$\delta a_i = 1 - a(t_i)$")

    ax5.scatter(range(ISImax), ISIs[:ISImax], s=10, fc="w", ec="k")
    ax5.scatter(range(ISImax), Tis_theo[1:ISImax + 1], marker="^", s=10, fc="w", ec="C1", zorder=1)
    ax5.axhline(T_infty, lw = 1, ls=":", c="C7")
    ax5.axhline(T0, lw = 1, ls=":", c="C7")
    ax5.text(3, T_infty-1, "$T_\infty$", va="top")
    ax5.text(7, T0, "$T_0$", va="bottom")
    ax5.set_xlabel("$i$")
    ax5.set_ylabel("$T_i$")

    difISIs = [T_infty - isi for isi in ISIs]

    Ts = [sum(ISIs[0:i]) for i in range(1, ISImax + 1)]
    ax2.scatter(Ts, ISIs[:ISImax], s=10, fc="w", ec="k")
    ax2.scatter([sum(Tis_theo[0:i]) for i in range(1, ISImax + 1)], Tis_theo[1:ISImax + 1], marker="^", s=10, fc="w", ec="C1", zorder=2)
    ax2.set_ylabel("$T_i(t)$")
    ax2.set_xlabel("$t$")
    ax2.axhline(T_infty, lw = 1, ls=":", c="C7")
    ax2.axhline(T0, lw = 1, ls=":", c="C7")
    ax2.text(100, T_infty-1, "$T_\infty$", va="top")
    ax2.text(250, T0, "$T_0$", va="bottom")
    ax2.arrow(600, T0, 0, T_infty - T0, fc="k", length_includes_head=True, head_width=20,
              head_length=1.0, lw=0.75, clip_on=False, zorder=5)
    ax2.arrow(600, T_infty, 0, -(T_infty - T0), fc="k", length_includes_head=True, head_width=20,
              head_length=1.0, lw=0.75, clip_on=False, zorder=5)
    ax2.text(550, (T0 + T_infty) / 2, r"$\Delta_a \tau_a$", ha="right")

    taueff = taua/(T_infty - taua*np.log(1-ampa))
    Ti_theo = []
    for i in range(ISImax):
        Ti_theo.append(T_infty - (T_infty - T0)*np.exp(-i/taueff))
    ax5.plot(range(ISImax), Ti_theo, c="C3")


    ax3.scatter(difISIs[0:ISImax], difISIs[1:ISImax + 1], s=10, fc="w", ec="k")
    ax3.scatter([Tis_theo[-1] - t for t in Tis_theo[1:ISImax]], [Tis_theo[-1] - t for t in Tis_theo[2:ISImax + 1]], marker="^", s=10, fc="w", ec="C1", zorder=2)
    ax3.set_xlabel(r"$T_\infty - T_i$")
    ax3.set_ylabel(r"$T_\infty - T_{i+1}$")



    popt, pcov = curve_fit(lin_func, difISIs[0:ISImax], difISIs[1:ISImax + 1], p0=(1))
    xs = np.linspace(0, 15, 100)
    ys = popt[0]*xs
    print("slope: ", popt[0], np.exp(-T_infty/taua)*(1-ampa))
    ax3.plot(xs, ys, ls="--", c="C7", zorder=3)
    ys_theo = [np.exp(-T_infty/taua)*(1-ampa)*x for x in xs]
    ax3.plot(xs, ys, c="C3", zorder=3)
    ratioISIs = [(T_infty - I2)/(T_infty - I1) for I1, I2 in zip(ISIs[0:ISImax], ISIs[1:ISImax+1])]
    ratioISIs_i_theory = []
    for i in range(ISImax):
        ratioISIs_i_theory.append(np.exp(-T_infty/taua)*(1 - ampa) + ampa*np.power(np.exp(-T_infty/taua), i))

    ax6.set_ylim([0, 1])
    ax6.scatter(range(ISImax), ratioISIs, s=10, fc="w", ec="k")
    ax6.scatter(range(ISImax-1), [(Tis_theo[-1] - I2) / (Tis_theo[-1] - I1) for I1, I2 in zip(Tis_theo[1:ISImax], Tis_theo[2:ISImax + 1])], marker="^", s=10, fc="w", ec="C1", zorder=2)
    ax6.axhline(np.exp(-T_infty / taua) * (1 - ampa), c="C3")
    ax6.set_ylabel(r"$\delta T_{i+1}/ \delta T_i$")
    ax6.set_xlabel(r"$i$")
    plt.savefig(home + "/Data/Calcium/Plots/13_map_ISIs.pdf", transparent=True)
    plt.show()


    print("T_\infty: ", T_infty)
    print("T0: ", T0)
    print(f"$T_\infty - T_0={T_infty - T0}$")
    print(f"$T_\infty/T_0={T_infty/T0}$")
    ratioT = (1 - ampa)*(1-np.exp(-T_infty/taua))/(1-(1-ampa)*np.exp(-T_infty/taua))
    print(rf"$T_\infty/T_0 \approx {1/ratioT}$")