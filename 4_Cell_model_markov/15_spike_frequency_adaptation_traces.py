import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from scipy.optimize import curve_fit

import functions as fc


def transient_func2(t, T0, T, tau):
    return T0*np.exp(-t/tau) + T*(1 - np.exp(-t/tau))


if __name__ == "__main__":
    tau = 3.00
    j = 0.05
    taua = 100
    ampa = 0.06
    home = os.path.expanduser("~")

    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    ISIs = np.loadtxt(home + folder + file_spike)

    RATEs = [1/isi for isi in ISIs[:100]]
    t_ISIs = []
    T = 0
    for isi in ISIs[:100]:
        T += isi
        t_ISIs.append(T)

    popt, pcov = curve_fit(transient_func2, t_ISIs, RATEs, p0=(10, 100, 200))

    fc.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    fc.remove_top_right_axis([ax1, ax2])

    t_plot = []
    a_plot = []
    for t, ca, jpuff, a in data:
        if t < max(t_ISIs):
            t_plot.append(t)
            a_plot.append(a)

    ax2.set_xlabel("$t$ / [s]")
    ax1.plot(t_plot, a_plot)

    ax2.set_xlabel("$t$ / [s]")
    ax2.set_ylabel(r"$f(t)$ / [1/s]")
    ax2.scatter(t_ISIs, RATEs, s=20, fc="w", ec="C0")

    xts = np.linspace(0, 1000)
    # yTs = popt[1]*(1 - np.exp(-(xts-popt[0])/popt[2]))
    yTs = popt[0] * np.exp(-xts / popt[2]) + popt[1] * (1 - np.exp(-xts / popt[2]))
    ax2.plot(xts, yTs, c="k", lw=1,
             label=f"$f_0$ = {popt[0]:.1e}" + "\n" + f"$f_\infty$ = {popt[1]:.1e}" + "\n" + rf"$\tau_a$ = {popt[2]:.1f}")
    legend = ax2.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 10})
    legend.get_frame().set_linewidth(0.5)

    plt.show()