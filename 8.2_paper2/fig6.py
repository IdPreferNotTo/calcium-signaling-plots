import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc


def transient_func(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))


if __name__ == "__main__":

    tau = np.logspace(-1, 2, 50)[33]
    jca = np.logspace(-3, 0, 50)[19]
    amp_a = 0.1
    tau_a = 1_000
    home = os.path.expanduser("~")
    folder = home + "/CLionProjects/PhD/calcium_spikes_langevin/out/Data_adap"
    file = f"/spike_times_langevin_ip1.00_taua{tau_a:.2e}_ampa{amp_a:.2e}_tau{tau:.2e}_j{jca:.2e}_N10_0.dat"

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(3, 1.5))
    gs = gridspec.GridSpec(1, 1)
    ax0 = fig.add_subplot(gs[0])
    st.remove_top_right_axis([ax0])

    ISIs = np.loadtxt(folder + file)
    ax0.scatter(range(len(ISIs)), ISIs)

    nr_ISIs = len(ISIs)
    index_ISIs = np.arange(nr_ISIs)

    popt, pcov = curve_fit(transient_func, index_ISIs, ISIs, p0=(50, 100, 2))
    i_ISI_fit = np.linspace(0, nr_ISIs)
    ISI_fit = popt[0] * np.exp(-i_ISI_fit / popt[2]) + popt[1] * (1 - np.exp(-i_ISI_fit / popt[2]))

    ax0.plot(i_ISI_fit, ISI_fit, c="k", zorder=5)
    print(popt[2])
    plt.show()