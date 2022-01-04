import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from scipy.optimize import curve_fit

import functions as fc
import styles as st

def transient_func2(t, T0, T, tau):
    return T0*np.exp(-t/tau) + T*(1 - np.exp(-t/tau))


if __name__ == "__main__":
    taua = 829
    ampa = 0.0309
    tau = 10.5
    j = 0.0146
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spikes = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    data = np.loadtxt(folder + file)
    ISIs = np.loadtxt(folder + file_spikes)

    RATEs = [1/isi for isi in ISIs[:100]]
    t_ISIs = []
    T = 0
    for isi in ISIs[:100]:
        t_ISIs.append(T)
        T += isi

    popt, pcov = curve_fit(transient_func2, t_ISIs, RATEs, p0=(10, 100, 200))
    f0 = popt[0]
    fi = popt[1]
    tauf = popt[2]
    perr = np.sqrt(np.diag(pcov))
    print(perr)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1:])
    st.remove_top_right_axis([ax1, ax2])

    t_plot = []
    a_plot = []
    for t, ca, jpuff, a in data:
        if t < max(t_ISIs):
            t_plot.append(t)
            a_plot.append(a)

    ax1.set_ylabel("$1-a(t)$")
    ax1.plot(t_plot, [1 - a for a in a_plot])
    ts = np.linspace(0, 1000, 5000)
    dt = ts[1] - ts[0]
    adaps = []
    bdaps = []
    cdaps = []
    b_ = 0
    r0 = np.loadtxt(home + "/Data/Calcium/data/r0_ip1.00_tau1.05e+01_j1.46e-02_N10.dat")
    bi = taua*ampa*fi/(1+taua*ampa*fi)
    b_theos = []
    b2_theos = []
    beta = 2*ampa*r0 + 1./taua
    alpha = ampa*r0
    b_root1 = beta/(2*alpha) + np.sqrt(np.power(beta/(2*alpha),2) -1)
    b_root2 = beta / (2 * alpha) - np.sqrt(np.power(beta / (2 * alpha), 2) - 1)
    print("r0: ", r0)
    for t in ts:
        bdaps.append(b_)
        b_ += ((-b_/taua) + np.power(1 - b_, 2)*ampa*r0)*dt
        b = b_root1*(1 - np.exp(-np.sqrt(4*taua*ampa*r0 +1)*t/taua))/(b_root1/b_root2 - np.exp(-np.sqrt(4*taua*ampa*r0 +1)*t/taua))
        b2 = b_root2*(1 - np.exp(-np.sqrt(4*taua*ampa*r0 +1)*t/taua))
        b_theos.append(b)
        b2_theos.append(b2)
    ax1.plot(ts, b_theos)
    ax1.plot(ts, b2_theos)
    ax1.axhline(1 - np.sqrt(4*taua*ampa*r0 +1)/(2*taua*ampa*r0) + 1/(2*taua*ampa*r0))
    t_eff =  taua/(1 + 2*taua*ampa*r0)
    print("t_eff: ", t_eff)

    ax2.set_xlabel("$t$ / [s]")
    ax2.set_ylabel(r"$f(t)$ / [1/s]")
    ax2.scatter(t_ISIs, RATEs, s=20, fc="w", ec="C0")

    xts = np.linspace(0, 1000)
    # yTs = popt[1]*(1 - np.exp(-(xts-popt[0])/popt[2]))
    yTs = popt[0] * np.exp(-xts / popt[2]) + popt[1] * (1 - np.exp(-xts / popt[2]))
    taua = tauf * 2* np.sqrt((r0/fi)*(r0 -fi)/fi)
    delta = (1/(fi*taua))*(r0 - fi)/fi
    print("fi/f0: ", fi/r0)
    print("1-b: ", (1-bi))
    print(taua, delta)
    ax2.plot(xts, yTs, c="k", lw=1,
             label=f"$f_0$ = {popt[0]:.1e} $\pm$ {perr[0]:.1e}" + "\n" + f"$f_\infty$ = {popt[1]:.1e} $\pm$ {perr[1]:.1e}" + "\n" + rf"$\tau_f$ = {popt[2]:.1f} $\pm$ {perr[2]:.1f}")
    legend = ax2.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 10})
    legend.get_frame().set_linewidth(0.5)

    ax1.set_xlim([0, 1000])
    ax2.set_xlim([0, 1000])
    plt.show()