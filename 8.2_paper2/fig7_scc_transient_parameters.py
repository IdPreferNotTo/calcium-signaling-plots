import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc


def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))


if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    axis = [ax1, ax2, ax3, ax4]
    st.remove_top_right_axis(axis)

    ax1.text(0.1, 0.95, "A$_{i}$", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.1, 0.95, "A$_{ii}$", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.1, 0.95, "B$_{i}$", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.1, 0.95, "B$_{ii}$", fontsize=13, transform=ax4.transAxes, va='top')
    for ax in axis:
        ax.set_ylim([-0.4, 0.1])
        ax.axhline(0, ls=":", c="C7")

    home = os.path.expanduser("~")
    folder = home + "/Data/calcium_spikes_langevin/Data_adap"
    folder_no_adap = home + "/Data/calcium_spikes_langevin/Data_no_adap"
    tau = np.logspace(-1, 2, 50)[33]
    jca = np.logspace(-3, 0, 50)[19]
    amp_as = np.logspace(-2, 0, 41)

    # VARY DELTA_A FIX TAU_A
    rho1s = []
    dTs = []
    iota_effs_s = []
    iota_effs_t = []
    a_inftys = []

    tau_as = np.logspace(1, 3, 41)
    tau_a = tau_as[35]
    print(tau_a)
    for amp_a in amp_as:
        file = f"/spike_times_langevin_ip1.00_taua{tau_a:.2e}_ampa{amp_a:.2e}_tau{tau:.2e}_j{jca:.2e}_N10_0.dat"
        file_no_adap = f"/spike_times_langevin_ip1.00_tau{tau:.2e}_j{jca:.2e}_N10_0.dat"
        ISIs = np.loadtxt(folder + file)
        ISIs_no_adap = np.loadtxt(folder_no_adap + file_no_adap)
        nr_ISIs = len(ISIs)
        index_ISIs = np.arange(nr_ISIs)

        popt, pcov = curve_fit(transient_func, index_ISIs, ISIs, p0=(50, 100, 2))
        i_ISI_fit = np.linspace(0, nr_ISIs)
        ISI_fit = popt[0] * np.exp(-i_ISI_fit / popt[2]) + popt[1] * (1 - np.exp(-i_ISI_fit / popt[2]))

        print(np.mean(ISIs[500:]), tau_a*np.log(1 - amp_a))
        iota_eff = tau_a /(np.mean(ISIs[500:]) - tau_a*np.log(1 - amp_a))
        iota_effs_t.append(iota_eff)
        iota_effs_s.append(popt[2])
        rho1 = fc.k_corr(ISIs[500:], ISIs[500:], 1) / fc.k_corr(ISIs[500:], ISIs[500:], 0)
        rho1s.append(rho1)
        dT = np.mean(ISIs[500:]) - np.mean(ISIs_no_adap)
        dTs.append(dT)
        a_infty = amp_a / (1. - (1 - amp_a)*np.exp(-popt[1]/tau_a))
        a_inftys.append(a_infty)

    cmap_viridis = plt.get_cmap("cividis", 31)


    ax1.set_xlabel(r"$\Delta T$")
    ax1.set_ylabel(r"$\rho_1$")
    im0 = ax1.scatter(dTs, rho1s, s=20, c=tau_as, cmap=cmap_viridis, norm=colors.LogNorm())

    ax2.set_xlim([0, 10])
    ax2.set_xlabel(r"$\iota_{\rm eff}$")
    ax2.set_ylabel("")
    ax2.plot(iota_effs_t, rho1s, c="k")
    cmap_viridis = plt.get_cmap("cividis", 31)
    im1 = ax2.scatter(iota_effs_s, rho1s, s=20, c=tau_as, cmap=cmap_viridis, norm=colors.LogNorm())
    cbar = fig.colorbar(im1, ax=ax2, ticks=[0.01, 0.1, 1])
    cbar.set_label(r"$\Delta_a$")

    # VARY TAU_A FIX DELTA_A

    rho1s = []
    dTs = []
    iota_effs_s = []
    iota_effs_t = []
    a_inftys = []

    amp_a = 0.1
    tau_as = np.logspace(1, 3, 41)
    for tau_a in tau_as:
        file = f"/spike_times_langevin_ip1.00_taua{tau_a:.2e}_ampa{amp_a:.2e}_tau{tau:.2e}_j{jca:.2e}_N10_0.dat"
        file_no_adap = f"/spike_times_langevin_ip1.00_tau{tau:.2e}_j{jca:.2e}_N10_0.dat"
        ISIs = np.loadtxt(folder + file)
        ISIs_no_adap = np.loadtxt(folder_no_adap + file_no_adap)
        nr_ISIs = len(ISIs)
        index_ISIs = np.arange(nr_ISIs)

        popt, pcov = curve_fit(transient_func, index_ISIs, ISIs, p0=(50, 100, 2))
        i_ISI_fit = np.linspace(0, nr_ISIs)
        ISI_fit = popt[0] * np.exp(-i_ISI_fit / popt[2]) + popt[1] * (1 - np.exp(-i_ISI_fit / popt[2]))

        print(np.mean(ISIs[500:]), tau_a*np.log(1 - amp_a))
        iota_eff = tau_a /(np.mean(ISIs[500:]) - tau_a*np.log(1 - amp_a))
        iota_effs_t.append(iota_eff)
        iota_effs_s.append(popt[2])
        rho1 = fc.k_corr(ISIs[500:], ISIs[500:], 1) / fc.k_corr(ISIs[500:], ISIs[500:], 0)
        rho1s.append(rho1)
        dT = np.mean(ISIs[500:]) - np.mean(ISIs_no_adap)
        dTs.append(dT)
        a_infty = amp_a / (1. - (1 - amp_a)*np.exp(-popt[1]/tau_a))
        a_inftys.append(a_infty)

    cmap_viridis = plt.get_cmap("cividis", 31)

    ax3.set_xlabel(r"$\Delta T$")
    ax3.set_ylabel(r"$\rho_1$")
    im2 = ax3.scatter(dTs, rho1s, s=20, c=tau_as, cmap=cmap_viridis, norm=colors.LogNorm())

    ax4.set_xlim([0, 5])
    ax4.set_xlabel(r"$\iota_{\rm eff}$")
    ax4.plot(iota_effs_t, rho1s, c="k")
    cmap_viridis = plt.get_cmap("cividis", 31)
    im3 = ax4.scatter(iota_effs_s, rho1s, s=20, c=tau_as, cmap=cmap_viridis, norm=colors.LogNorm())
    cbar = fig.colorbar(im3, ax=ax4, ticks=[10, 100, 1000])
    cbar.set_label(r"$\tau_a$")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig6.png", transparent=True)
    plt.show()