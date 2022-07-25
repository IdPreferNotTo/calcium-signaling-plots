import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings
import os

import styles as st


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
    taus = [4.47, 1.26]
    js = [0.0178, 0.0562]
    IP3s = np.linspace(0.02, 2, 100)
    IP3s_density = [[IP3s[49], IP3s[53], IP3s[64]], [IP3s[49], IP3s[50], IP3s[52]]]

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
        st.remove_top_right_axis([ax1, ax2, ax3])
        folder_markov_ip3 = home + "/Data/calcium_spikes_markov/Data_no_adap_ip3/"
        folder_langevin_ip3 = home + "/Data/calcium_spikes_langevin_strat/Data_no_adap_ip3/"

        print("Load mean, cv over ip3...")
        rates_markov = []
        cvs_markov = []
        file_theory_r0 = home + f"/Data/calcium_spikes_theory/ca_langevin_strat_r0_tau{tau:.2e}_j{j:.2e}_over_ip3.dat"
        data_theory_r0 = np.loadtxt(file_theory_r0)
        ip3s_theory, r0_theory = np.transpose(data_theory_r0)
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

        ax1.set_xlabel("$s$")
        ax1.set_ylabel(r"$r_0$ / s$^{-1}$")
        ax1.set_xlim([0.5, 1.5])
        if k == 0:
            ax1.set_ylim([0, 0.11])
        else:
            ax1.set_ylim([0., 0.22])
        ax1.plot(ip3s_theory, r0_theory, color=color.palette[5], label="Theory")

        cmap_viridis = plt.get_cmap("cividis", 11)
        im1 = ax1.scatter(IP3s[::2], rates_markov[::2], s=20, c=cvs_markov[::2], cmap=cmap_viridis, vmin=-0.05,
                          vmax=1.05, label="Sim.")
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
        cbar.ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        cbar.set_label(r"$CV_T$", loc="center")

        if k == 0:
            ax1.legend(fancybox=False, fontsize=10, loc=1)

        isis_markovs = []
        print("Load ISI densities...")
        if k == 0:
            plt_color = color.palette[1]
        else:
            plt_color = color.orange[1]
        for n, IP3 in enumerate(IP3s_density[k]):
            p_isis = []
            file_spikes_markov = f"spike_times_markov_ip{IP3:.2f}_tau{tau:.2e}_j{j:.2e}_K10_5.dat"
            file_fano_markov = f"/Data/calcium_spikes_theory/Fano_factor_ip{IP3:.2f}_tau{tau:.2e}_j{j:.2e}.dat"
            warnings.simplefilter("ignore")
            isis_markov = np.loadtxt(folder_markov_ip3 + file_spikes_markov)
            isis_markovs.append(isis_markov)
            mean_isi = np.mean(isis_markov)
            std_isi = np.std(isis_markov)

            cv_isi = std_isi / mean_isi
            ts_inv_gau = np.linspace(0.1, 100, 1001)
            inv_gaus1 = []
            for t in ts_inv_gau:
                p = np.sqrt(mean_isi / (2 * np.pi * np.power(cv_isi, 2) * (t ** 3))) * np.exp(
                    -(t - mean_isi) ** 2 / (2 * mean_isi * np.power(cv_isi, 2) * t))
                inv_gaus1.append(p)
            if k == 0:
                bins = 20
            else:
                bins = 50
            ax2.hist(isis_markov, bins=bins, color=plt_color, density=True, alpha=0.2 + 0.3 * (2 - n),
                     label=f"$s = {IP3:.2f}$")
            ax2.plot(ts_inv_gau, inv_gaus1, ls="--", color=color.palette[5], alpha=0.4 + 0.3 * (2 - n))

            data_fano = np.loadtxt(home + file_fano_markov)
            deltaT, mean, var, fano = np.transpose(data_fano)
            ax3.plot(deltaT, fano, color=plt_color, alpha=0.2 + 0.3 * (2 - n))
            ax3.axvline(mean_isi, color=plt_color, ls="--", alpha=0.2 + 0.3 * (2 - n))
            if n == k and k == 0:
                ax3.axvline(mean_isi, color=plt_color, ls="--", alpha=0.2 + 0.3 * (2 - n), label=r"$\langle T \rangle$")
                ax3.legend(fancybox=False, fontsize=9, loc=1)

            ax2.set_ylim([0, 1.3 * max(inv_gaus1)])
        ax2.legend(fancybox=False, fontsize=9, loc=1)
        ax2.set_ylabel(r"$p_{\rm ISI}(T)$")
        ax2.set_xlabel("$T$ / s")

        ax3.set_xlim([1, 1000])
        ax2.set_xlim([0, 100])

        if k == 1:
            ax2.set_ylim([0, 0.11])
        ax3.set_xscale("log")
        ax3.set_xlabel("$t$ / s")
        ax3.set_ylabel("$F(t)$")
        ax3.set_ylim([0, 1.3])

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig7b.pdf")
    plt.show()
