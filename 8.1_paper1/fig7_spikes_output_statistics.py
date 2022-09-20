import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import warnings
import os

import styles as st
import functions as fc


def get_inverse_gaussian(ts, mean, cv):
    inv_gaus = []
    for t in ts:
        p = np.sqrt(mean / (2 * np.pi * np.power(cv, 2) * (t ** 3))) * np.exp(
            -(t - mean) ** 2 / (2 * mean * np.power(cv, 2) * t))
        inv_gaus.append(p)
    return inv_gaus


if __name__ == "__main__":
    # Parameters
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(9, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=3)
    color = st.Colors

    home = os.path.expanduser("~")
    taus = [5.62, 1.78]
    js = [0.0126, 0.0355]
    IP3s = np.linspace(0.50, 1.50, 101)
    IP3s = IP3s[:-1]

    for k, (tau, j) in enumerate(zip(taus, js)):
        ax1 = fig.add_subplot(gs[k, 0])
        ax2 = fig.add_subplot(gs[k, 1])
        ax3 = fig.add_subplot(gs[k, 2])
        if k == 0:
            label = "A"
        else:
            label = "B"
        ax1.text(0.05, 0.95, label + r"$_{\rm i}$", fontsize=13, transform=ax1.transAxes, va='top')
        ax2.text(0.05, 0.95, label + r"$_{\rm ii}$", fontsize=13, transform=ax2.transAxes, va='top')
        ax3.text(0.05, 0.95, label + r"$_{\rm iii}$", fontsize=13, transform=ax3.transAxes, va='top')
        st.remove_top_right_axis([ax1, ax3, ax2])
        folder_markov_ip3 = home + "/Data/calcium_spikes_markov/Data_no_adap_ip3/"
        folder_langevin_ip3 = home + "/Data/calcium_spikes_langevin_strat/Data_no_adap_ip3/"

        print("Load mean, cv over ip3...")
        rates_markov = []
        cvs_markov = []
        file_theory_r0 = home + f"/Data/calcium_spikes_theory/ca_langevin_r0_tau{tau:.2e}_j{j:.2e}_over_ip3.dat"
        data_theory_r0 = np.loadtxt(file_theory_r0)
        ip3s_theory, r0_theory = np.transpose(data_theory_r0)
        r0_theory2 = []
        for ip3 in ip3s_theory:
            file_langevin_r0 = folder_langevin_ip3 + f"/spike_times_langevin_ip{ip3:.2f}_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
            isis_langevin = np.loadtxt(file_langevin_r0)
            mean_isis_langevin = np.mean(isis_langevin)
            r0_theory2.append(1/mean_isis_langevin)
        print(ip3s_theory)
        for IP3 in IP3s:
            file_markov = f"spike_times_markov_ip{IP3:.2f}_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
            warnings.simplefilter("ignore")
            isis_markov = np.loadtxt(folder_markov_ip3 + file_markov)
            if len(isis_markov) == 0:
                rates_markov.append(0)
                cvs_markov.append(1)
            else:
                mean_isi_markov = np.mean(isis_markov)
                std_isi_markov = np.std(isis_markov)
                rates_markov.append(1. / mean_isi_markov)
                cvs_markov.append(std_isi_markov / mean_isi_markov)

        if k == 0:
            plt_color = color.palette[0]
        else:
            plt_color = color.palette[2]
        if k == 0:
            arrow = mpatches.FancyArrow(1, rates_markov[51] + 0.04, 0, -0.03, width=0.005, head_length=0.025/2,
                                        head_width= 0.03, length_includes_head=True, color=plt_color)
        else:
            arrow = mpatches.FancyArrow(1, rates_markov[51] + 0.08, 0, -0.06, width=0.005, head_length=0.05/2,
                                        head_width= 0.03, length_includes_head=True, color=plt_color)
        ax1.add_patch(arrow)

        ax1.set_xlabel("$s$")
        ax1.set_ylabel(r"$r_0$ / s$^{-1}$")
        ax1.set_xlim([0.5, 1.5])
        if k == 0:
            ax1.set_ylim([0, 0.11])
        else:
            ax1.set_ylim([0., 0.22])
        ax1.plot(ip3s_theory, r0_theory2, color=color.palette[5], label="Theory")

        cmap_viridis = plt.get_cmap("cividis", 11)
        im1 = ax1.scatter(IP3s[1::2], rates_markov[1::2], s=20, c=cvs_markov[1::2], cmap=cmap_viridis, vmin=-0.05,
                          vmax=1.05, label="Sim.")
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
        cbar.ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        cbar.set_label(r"$CV_T$", loc="center")

        if k == 0:
            ax1.legend(fancybox=False, fontsize=10, loc=1, framealpha=1.0)

        isis_markovs = []
        print("Load ISI densities...")
        p_isis = []
        file_spikes_markov = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
        file_fano_markov = f"/Data/calcium_spikes_theory/Fano_factor_ip1.00_tau{tau:.2e}_j{j:.2e}.dat"
        warnings.simplefilter("ignore")
        isis_markov = np.loadtxt(folder_markov_ip3 + file_spikes_markov)
        isis_markovs.append(isis_markov)
        mean_isi = np.mean(isis_markov)
        std_isi = np.std(isis_markov)
        cv_isi = std_isi / mean_isi
        print(cv_isi)

        fs = np.logspace(-3, 0, 100)
        spectrum_data = fc.power_spectrum_isis(fs, isis_markov, Tmax=1000)
        spectrum_theory = []
        for f in fs:
            cv2 = np.power(cv_isi, 2)
            r0 = np.power(mean_isi, -1)
            arg = (1 - np.sqrt(1- 1j*4*np.pi*f*mean_isi*cv2))/cv2
            sf = r0*(1 - np.power(np.absolute(np.exp(arg)),2))/(np.power(np.absolute(1-np.exp(arg)), 2))
            spectrum_theory.append(sf)

        ax3.set_xlabel("$f$ / s$^{-1}$")
        ax3.set_ylabel("$S(f)$")
        ax3.axhline(1/mean_isi, ls=":", color="C7")
        ax3.axhline((cv_isi)**2/mean_isi, ls=":", color="C7")
        ax3.set_xlim([10**(-3), 1])
        ax3.set_ylim([0.0004, 0.5])
        ax3.set_yscale("log")
        ax3.set_xscale("log")
        ax3.plot(fs, spectrum_data, color=plt_color)
        ax3.plot(fs, spectrum_theory, color=color.palette[5], lw=1, ls="--")
        ax3.text(0.003, 1.6/mean_isi, "$r_0$")
        if k==0:
            ax3.axvline(1/mean_isi, ls=":", color="C7", zorder=1)

        if k == 0:
            ax3.text(0.1, 1.9*cv_isi**2/mean_isi, "$r_0 CV_T^2$")
        else:
            ax3.text(0.1, 0.35 * cv_isi ** 2 / mean_isi, "$r_0 CV_T^2$")

        data_fano = np.loadtxt(home + file_fano_markov)
        deltaT, mean, var, fano = np.transpose(data_fano)
        ts = np.linspace(0, mean_isi, 100)
        fano_det = []
        for t in ts:
            fano_det.append(1 - t/mean_isi)
        #ax2.plot(ts, fano_det, c="k", ls="--", zorder=3)
        ax2.plot(deltaT, fano, color=plt_color, zorder=2)
        ax2.axhline(cv_isi**2, ls=":", color="C7", zorder=1)
        ax2.axvline(mean_isi, ls=":", color="C7", zorder=1)

        ax2.set_xlim([1, 1000])
        ax2.set_xscale("log")
        ax2.set_xlabel("$t$ / s")
        ax2.set_ylabel("$F(t)$")
        ax2.set_ylim([0, 1.3])
        ax2.text(70, 1.1, r"$\langle T \rangle$")
        ax2.text(3., 0.1 + cv_isi**2, "$CV_T^2$")

plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig7.pdf")
plt.show()
