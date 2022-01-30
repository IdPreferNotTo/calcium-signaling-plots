import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import os
import styles as st


def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))


def get_isis_from_data(data, nr_cell):
    ts = [x[0] for x in data]
    cas = [x[nr_cell] for x in data]
    spiking: bool = False
    t_tmp: float = 0
    spike_times = []
    if nr_cell == 2:
        for t, ca in zip(ts, cas):
            if t < 3700:
                if ca > 0.35 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.35 and spiking:
                    spiking = False
    else:
        for t, ca in zip(ts, cas):
            if t < 3700:
                if ca > 0.4 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.4 and spiking:
                    spiking = False
    ISIs = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]
    t_ISIs = []
    T = 0
    for isi in ISIs:
        T += isi
        t_ISIs.append(T)
    return ISIs


if __name__ == "__main__":
    home = os.path.expanduser("~")
    file = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
    data = np.loadtxt(file)
    nr_cell = 13 #2, 12
    isis_exp = get_isis_from_data(data, nr_cell)
    nr_isis_exp = len(isis_exp)
    i_isis_exp = np.arange(nr_isis_exp)
    popt_exp, pcov = curve_fit(transient_func, i_isis_exp, isis_exp, p0=(100, 150, 2))
    print(popt_exp)
    isis_exp_fit = popt_exp[0]*np.exp(-i_isis_exp/popt_exp[2]) + popt_exp[1]*(1 - np.exp(-i_isis_exp/popt_exp[2]))

    taua = 829
    ampa = 0.0309
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap"
    file_isi = f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    isis_sim = np.loadtxt(folder + file_isi)
    popt_sim, pcov = curve_fit(transient_func, i_isis_exp, isis_sim[:nr_isis_exp], p0=(100, 150, 2))
    print(popt_sim)
    isis_sim_fit = popt_sim[0] * np.exp(-i_isis_exp / popt_sim[2]) + popt_sim[1] * (1 - np.exp(-i_isis_exp / popt_sim[2]))

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3])

    # Data
    ax0.set_xlim([0, nr_isis_exp])
    ax0.set_ylim([0, 1.5*popt_exp[1]])
    ax0.set_xlabel("$i$")
    ax0.set_ylabel("$T_i$")
    ax0.scatter(range(nr_isis_exp), isis_exp, ec=st.colors[4], fc="w", s=20, zorder=3)
    ax0.plot(i_isis_exp, isis_exp_fit, lw=1, c="k", zorder=2, label = "fit")
    ax0.axhline(popt_exp[0], ls=":", c="C7")
    ax0.axhline(popt_exp[1], ls=":", c="C7")
    ax0.text(len(isis_exp)-10, popt_exp[0]*1.1, "$T_0$", ha="center")
    ax0.text(5, popt_exp[1]*1.1, "$T_\infty$", ha="center")

    ax1.set_xlim([-35, 110])
    ax1.set_ylim([-35, 110])
    ax1.set_xlabel("$T_\infty - T_{i}$")
    ax1.set_ylabel("$T_\infty - T_{i+1}$")
    ax1.scatter([popt_exp[1] - I for I in isis_exp[:-1]], [popt_exp[1] - I for I in isis_exp[1:]], ec=st.colors[4], fc="w", s=20, zorder=3)
    ax1.axhline(0, ls=":", c="C7")
    ax1.axvline(0, ls=":", c="C7")
    m_sim = np.exp(-1/popt_exp[2])
    xs = np.linspace(-25, 75, 200)
    ys = [m_sim*x for x in xs]
    ax1.plot(xs, ys, lw=1, c="k", zorder=2)
    ax1.text(5, 90, r"$\alpha = e^{-1/\iota_{\rm eff}}$", ha="left", va="center", fontsize = 15)

    # Simulation
    ax2.set_xlim([0, nr_isis_exp])
    ax2.set_ylim([0, 1.5*popt_exp[1]])
    ax2.set_xlabel("$i$")
    ax2.set_ylabel("$T_i$")
    ax2.scatter(range(nr_isis_exp), isis_sim[:nr_isis_exp], ec=st.colors[1], fc="w", s=20, zorder=3)
    ax2.plot(i_isis_exp, isis_sim_fit, lw=1, c="k", zorder=2)
    ax2.axhline(popt_sim[0], ls=":", c="C7")
    T0 = popt_sim[0]
    ax2.axhline(popt_sim[1], ls=":", c="C7")
    ax2.text(len(isis_exp)-10, popt_sim[0]*1.1, "$T_0$", ha="center")
    ax2.text(5, popt_sim[1]*1.1, "$T_\infty$", ha="center")

    ax3.set_xlim([-35, 110])
    ax3.set_ylim([-35, 110])
    ax3.set_xlabel("$T_\infty - T_{i}$")
    ax3.set_ylabel("$T_\infty - T_{i+1}$")
    ax3.axhline(0, ls=":", c="C7")
    ax3.axvline(0, ls=":", c="C7")
    m_sim = np.exp(-1/popt_sim[2])
    xs = np.linspace(-25, 75, 200)
    ys = [m_sim*x for x in xs]
    m_theo = np.exp(-popt_exp[1]/taua)*(1.-ampa)
    ys2 = [m_theo * x for x in xs]
    print(m_sim, m_theo)
    ax3.scatter([popt_sim[1] - I for I in isis_sim[0:nr_isis_exp]], [popt_sim[1] - I for I in isis_sim[1:nr_isis_exp+1]], ec=st.colors[1], fc="w", s=20, zorder=3)
    ax3.plot(xs, ys, lw=1, c="k", zorder=2)
    ax3.plot(xs, ys2, lw=1, ls=":", c="k", zorder=2)
    ax3.text(5, 90, r"$\alpha = e^{-T_\infty/\tau_a}(1 - \Delta_a)$", ha="left", va="center", fontsize = 15)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig5.pdf", transparent=True)
    plt.show()