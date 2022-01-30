import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st


def inverse_gaussian(T, CV):
    ps = []
    ts = np.linspace(0, 2*T, 500)
    for t in ts:
        p = np.sqrt(T / (2 * np.pi * (CV**2) * (t ** 3))) * np.exp(-(t - T) ** 2 / (2 * T * (CV**2) * t))
        ps.append(p)
    return ts, ps

if __name__ == "__main__":
    taua1 = 829
    ampa1 = 0
    taua2 = 829
    ampa2 = np.logspace(-2, 1, 50)[10]
    taua3 = 829
    ampa3 = np.logspace(-2, 1, 50)[15]

    home = os.path.expanduser("~")
    folder1 = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_no_adap"
    file_ca_1 = f"/ca_markov_ip1.00_tau1.05e+01_j1.46e-02_N10_0.dat"
    file_isi_1 = f"/spike_times_markov_ip1.00_tau1.05e+01_j1.46e-02_N10_0.dat"
    data_ca_1 = np.loadtxt(folder1 + file_ca_1)
    data_ca_1 = [x for x in data_ca_1 if x[1]!=1]
    ts_1, cas_1, jpuffs_1, adaps_1 = np.transpose(data_ca_1)
    isis_1 = np.loadtxt(folder1 + file_isi_1)

    folder2 = home + "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap"
    file_ca_2 = f"/ca_markov_ip1.00_taua{taua2:.2e}_ampa{ampa2:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    file_isi_2 = f"/spike_times_markov_ip1.00_taua{taua2:.2e}_ampa{ampa2:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    data_ca_2 = np.loadtxt(folder2 + file_ca_2)
    data_ca_2 = [x for x in data_ca_2 if x[1] != 1]
    ts_2, cas_2, jpuffs_2, adaps_2 = np.transpose(data_ca_2)
    isis_2 = np.loadtxt(folder2 + file_isi_2)

    file_ca_3 = f"/ca_markov_ip1.00_taua{taua3:.2e}_ampa{ampa3:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    file_isi_3 = f"/spike_times_markov_ip1.00_taua{taua3:.2e}_ampa{ampa3:.2e}_tau1.05e+01_j1.46e-02_N10_0.dat"
    data_ca_3 = np.loadtxt(folder2 + file_ca_3)
    data_ca_3 = [x for x in data_ca_3 if x[1] != 1]
    ts_3, cas_3, jpuffs_3, adaps_3 = np.transpose(data_ca_3)
    isis_3 = np.loadtxt(folder2 + file_isi_3)


    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3*2.25, 2*2))
    gs = gridspec.GridSpec(2, 3)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[0, 2])
    ax5 = fig.add_subplot(gs[1, 2])
    st.remove_top_right_axis([ax0, ax1, ax2, ax3, ax4, ax5])

    ax0.text(0.1, 0.95, "A$_{i}$", fontsize=13, transform=ax0.transAxes, va='top')
    ax1.text(0.1, 0.95, "A$_{ii}$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, "B$_{i}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.1, 0.95, "B$_{ii}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.1, 0.95, "C$_{i}$", fontsize=13, transform=ax4.transAxes, va='top')
    ax5.text(0.1, 0.95, "C$_{ii}$", fontsize=13, transform=ax5.transAxes, va='top')


    ax0.set_ylim([0, 3.5])
    ax0.set_xlabel("$c_i$")
    ax0.set_ylabel("$P_0(c_i)$")
    ax0.hist(cas_1, bins=20, color=st.colors[1], density=True, alpha=0.6)
    file_p0_theory = home + f"/Data/Calcium/data/p0_ip1.00_tau1.05e+01_j1.46e-02_N10.dat"
    data = np.loadtxt(file_p0_theory)
    ca_theory, p0_theory = np.transpose(data)
    ax0.plot(ca_theory[::100], p0_theory[::100], color=st.colors[0])


    mean_isis_1 = np.mean(isis_1)
    std_isis_1 = np.std(isis_1)
    cv_isis_1 = std_isis_1/mean_isis_1
    ts_inv_gaussian, ps_inv_gaussian = inverse_gaussian(mean_isis_1, cv_isis_1)
    ax1.set_xlim(0, 2*mean_isis_1)
    ax1.set_xticks([0, mean_isis_1, 2*mean_isis_1])
    ax1.set_xticklabels(["$0$", r"$\langle T \rangle$", r"2$\langle T \rangle$"])
    ax1.set_xlabel("$T_i$")
    ax1.set_ylabel("$P_0(T_i)$")
    ax1.hist(isis_1, bins=20, color=st.colors[1], density=True, alpha=0.6)
    ax1.plot(ts_inv_gaussian, ps_inv_gaussian, color="k",  ls="--")

    ax2.set_ylim([0, 3.5])
    ax2.set_xlabel("$c_i$")
    ax2.set_ylabel("$P_0(c_i)$")
    ax2.hist(cas_2, bins=20, color=st.colors[1], density=True, alpha=0.6)

    mean_isis_2 = np.mean(isis_2)
    std_isis_2 = np.std(isis_2)
    cv_isis_2 = std_isis_2/mean_isis_2
    ts_inv_gaussian, ps_inv_gaussian = inverse_gaussian(mean_isis_2, cv_isis_2)
    ax3.set_xlim(0, 2*mean_isis_2)
    ax3.set_xticks([0, mean_isis_2, 2*mean_isis_2])
    ax3.set_xticklabels(["$0$", r"$\langle T \rangle$", r"2$\langle T \rangle$"])
    ax3.set_xlabel("$T_i$")
    ax3.set_ylabel("$P_0(T_i)$")
    ax3.hist(isis_2, bins=20, color=st.colors[1], density=True, alpha=0.6)
    ax3.plot(ts_inv_gaussian, ps_inv_gaussian, color="k", ls="--")

    ax4.set_ylim([0, 3.5])
    ax4.set_xlabel("$c_i$")
    ax4.set_ylabel("$P_0(c_i)$")
    ax4.hist(cas_3, bins=20, color=st.colors[1], density=True, alpha=0.6)

    mean_isis_3 = np.mean(isis_3)
    std_isis_3 = np.std(isis_3)
    cv_isis_3 = std_isis_3/mean_isis_3
    ts_inv_gaussian, ps_inv_gaussian = inverse_gaussian(mean_isis_3, cv_isis_3)
    ax5.set_xlim(0, 2*mean_isis_3)
    ax5.set_xticks([0, mean_isis_3, 2*mean_isis_3])
    ax5.set_xticklabels(["$0$", r"$\langle T \rangle$", r"2$\langle T \rangle$"])
    ax5.set_xlabel("$T_i$")
    ax5.set_ylabel("$P_0(T_i)$")
    ax5.hist(isis_3, bins=20, color=st.colors[1], density=True, alpha=0.6)
    ax5.plot(ts_inv_gaussian, ps_inv_gaussian, color="k", ls="--")

    print(ampa1, ampa2, ampa3)
    print(mean_isis_1, mean_isis_2, mean_isis_3)
    print(cv_isis_1, cv_isis_2, cv_isis_3)

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig2.png", transparent=True)
    plt.show()