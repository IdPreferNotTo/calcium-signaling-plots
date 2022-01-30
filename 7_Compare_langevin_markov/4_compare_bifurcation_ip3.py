import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    ip3s = np.linspace(0, 5, 100)
    rates_langevin = []
    CVs_langevin = []
    rates_markov= []
    CVs_markov = []
    n = 4
    tau = 5.96
    jca = 0.0596
    for ip3 in ip3s:
        folder_langevin: str = "/CLionProjects/PhD/calcium_spikes_langevin/out/"
        folder_markov: str = "/CLionProjects/PhD/calcium_spikes_markov/out/"
        file_langevin: str = "spike_times_langevin_ip{:.2f}_taua1.00e+02_ampa2.00e-01_tau{:.2e}_j{:.2e}_N10_{:d}.dat".format(ip3, tau, jca, 0)
        file_markov: str =  "spike_times_markov_ip{:.2f}_taua1.00e+02_ampa2.00e-01_tau{:.2e}_j{:.2e}_N10_{:d}.dat".format(ip3, tau, jca, 0)

        data_langevin = np.loadtxt(home + folder_langevin + file_langevin)
        data_markov = np.loadtxt(home + folder_markov + file_markov)
        if len(data_langevin) == 0:
            rate_langevin = 0
            CV_langevin = np.nan
        else:
            T = np.mean(data_langevin)
            std = np.std(data_langevin)
            rate_langevin = 1/T
            CV_langevin = std/T
        rates_langevin.append(rate_langevin)
        CVs_langevin.append(CV_langevin)

        if len(data_markov) == 0:
            rate_markov = 0
            CV_markov = np.nan
        else:
            T = np.mean(data_markov)
            std = np.std(data_markov)
            rate_markov = 1/T
            CV_markov = std/T
        rates_markov.append(rate_markov)
        CVs_markov.append(CV_markov)


    fig, ax = plt.subplots(1, 1)
    ax.plot(ip3s, rates_markov, c="k", label="Spike rate markov")
    ax.plot(ip3s, CVs_markov, c="k", ls=":", label="Spike CV markov")

    ax.plot(ip3s, rates_langevin, c="C7", label="Spike rate langevin")
    ax.plot(ip3s, CVs_langevin, c="C7", ls=":", label="Spike CV langevin")

    ax.set_xlabel("IP$_3$")
    ax.set_ylabel("Spike rate $r$")
    ax.set_ylim([-0.01, 0.2])
    ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
    ax.set_xlim([0, 5])
    secax= ax.secondary_yaxis(
        'right')
    secax.set_ylabel(r"Coefficient of variation $C_V$", rotation=-90)
    secax.set_yticks([0, 0.05, 0.1, 0.15, 0.2])

    ax.legend(prop={"size": 7}, loc=1, ncol=1, framealpha=1., edgecolor="k")
    leg = ax.get_legend()
    leg.get_frame().set_linewidth(0.5)


    axins = inset_axes(ax, width="40%", height="30%", loc='lower left',
                      bbox_to_anchor=(0.5, 0.6, 1, 1), bbox_transform=ax.transAxes)
    puff_rates = []
    for ip3 in ip3s:
        r_opn = n*0.13*np.power(ip3, 3)*(2./(1. + np.power(ip3, 3)))
        r_ref = n*1.3*np.power(ip3, 3)*(2./(1. + np.power(ip3, 3)))
        T = 1./r_opn + 1./r_ref
        puff_rates.append(1/T)
    axins.plot(ip3s, puff_rates, c="k")
    axins.set_xlim([0, 5])
    axins.set_xlabel("IP$_3$")
    axins.set_ylabel(r"Puff rate $r_{\rm puff}$")
    plt.show()