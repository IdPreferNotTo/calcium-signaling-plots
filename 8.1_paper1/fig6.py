import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import functions as fc
import styles as st

def get_inverse_gaussian(ts, mean, cv):
    inv_gaus = []
    for t in ts:
        p = np.sqrt(mean / (2 * np.pi * np.power(cv,2) * (t ** 3))) * np.exp(
            -(t - mean) ** 2 / (2 * mean * np.power(cv,2) * t))
        inv_gaus.append(p)
    return inv_gaus

if __name__ == "__main__":
    # Parameters
    P_ca = []
    P_isi = []
    P_ca_theory = []

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4.5, 3))
    grids = gridspec.GridSpec(1, 1)
    ax0 = fig.add_subplot(grids[0])
    axin1 = ax0.inset_axes([0.15, 0.65, 0.35, 0.35])
    axin2 = ax0.inset_axes([0.65, 0.15, 0.35, 0.35])
    axis = [ax0, axin1, axin2]
    st.remove_top_right_axis(axis)
    home = os.path.expanduser("~")

    ax0.set_xlabel("IP$_3$")
    ax0.set_xlim([0.0, 2.0])
    ax0.set_ylabel(r"$1/\langle T \rangle$ / s$^{-1}$")
    IP3s = np.linspace(0.02, 2, 100)
    meanISIs = []
    rates = []
    folder_ip3 = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_no_adap_ip3/"
    for IP3 in IP3s:
       file = f"spike_times_markov_ip{IP3:.2f}_tau1.05e+01_j1.10e-02_N10_0.dat"
       ISIs = np.loadtxt(folder_ip3 + file)
       if len(ISIs) == 0:
           rates.append(0)
       else:
           meanISI = np.mean(ISIs)
           meanISIs.append(meanISI)
           rates.append(1/meanISI)
    IP3s_theory = np.linspace(0.5, 2, 16)
    rates_theory = []
    for IP3_theory in IP3s_theory:
       print(f"{IP3_theory:.2f}")
       file_theory = home + f"/Data/Calcium/data/r0_ip{IP3_theory:.2f}_tau1.05e+01_j1.10e-02_N10.dat"
       rate = np.loadtxt(file_theory)
       rates_theory.append(rate)

    ax0.plot(IP3s, rates, lw=1, color=st.colors[1])
    ax0.plot(IP3s_theory, rates_theory, lw=1, color=st.colors[2])

    file = f"spike_times_markov_ip1.00_tau1.05e+01_j1.10e-02_N10_0.dat"
    ISIs_ip1 = np.loadtxt(folder_ip3 + file)
    mean_ISI_ip1 = np.mean(ISIs_ip1)
    cv_ISI_ip1 = np.std(ISIs_ip1)/mean_ISI_ip1
    axin1.hist(ISIs_ip1, bins=20, color=st.colors[1], density=True, alpha=0.6)
    axin1.set_xticks([0, mean_ISI_ip1, 2*mean_ISI_ip1])
    axin1.set_xticklabels(["$0$", r"$\langle T \rangle$", r"$2\langle T \rangle$"])
    axin1.set_ylabel("$P(T)$")
    axin1.set_yticks([])
    ts_inv_gau = np.linspace(0.1, 2*mean_ISI_ip1, 1001)
    inv_gaus = []
    for t in ts_inv_gau:
        p = np.sqrt(mean_ISI_ip1 / (2 * np.pi * np.power(cv_ISI_ip1, 2) * (t ** 3))) * np.exp(
            -(t - mean_ISI_ip1) ** 2 / (2 * mean_ISI_ip1 * np.power(cv_ISI_ip1, 2) * t))
        inv_gaus.append(p)
    axin1.plot(ts_inv_gau, inv_gaus, color="k", lw=1, label="Inv.\ Gaussian")
    axin1.text(0.05, 0.8, "IP$_3$=1", transform=axin1.transAxes)

    file = f"spike_times_markov_ip2.00_tau1.05e+01_j1.10e-02_N10_0.dat"
    ISIs_ip2 = np.loadtxt(folder_ip3 + file)
    mean_ISI_ip2 = np.mean(ISIs_ip2)
    cv_ISI_ip2 = np.std(ISIs_ip2)/mean_ISI_ip2
    axin2.hist(ISIs_ip2, bins=20, color=st.colors[1], density=True, alpha=0.6)
    axin2.set_xticks([0, mean_ISI_ip2, 2 * mean_ISI_ip2])
    axin2.set_xticklabels(["$0$", r"$\langle T \rangle$", r"$2\langle T \rangle$"])
    axin2.set_ylabel(r"$P(T)$")
    axin2.set_yticks([])
    ts_inv_gau = np.linspace(0.1, 2*mean_ISI_ip2, 1001)
    inv_gaus = []
    for t in ts_inv_gau:
        p = np.sqrt(mean_ISI_ip2 / (2 * np.pi * np.power(cv_ISI_ip2, 2) * (t ** 3))) * np.exp(
            -(t - mean_ISI_ip2) ** 2 / (2 * mean_ISI_ip2 * np.power(cv_ISI_ip2, 2) * t))
        inv_gaus.append(p)
    axin2.plot(ts_inv_gau, inv_gaus, color="k", lw=1, label="Inv.\ Gaussian")
    axin2.text(0.05, 0.8, "IP$_3$=2", transform=axin2.transAxes)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig6.pdf",transparent=True)
    plt.show()
