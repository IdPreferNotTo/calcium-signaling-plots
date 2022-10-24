import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

import functions as fc
import styles as st
import default_parameters as df


if __name__ == "__main__":
    home = os.path.expanduser("~")
    # Parameters
    tau = 1
    p = 0.06
    M = 4
    N = 5
    K = 10
    print(tau, p)
    # Load Data
    data_isi = df.load_spike_times_markov(tau, p, cer=False)

    mean_isi = np.mean(data_isi)
    std_isi = np.std(data_isi)

    ts = np.linspace(0, 2*mean_isi, 100)
    density_inverse_gaus = fc.inverse_gaussian_dist(ts, mean_isi, (std_isi / mean_isi) ** 2)
    density_gamma = fc.inverse_gaussian_dist(ts, mean_isi, (std_isi / mean_isi) ** 2)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    st.remove_top_right_axis([ax])


    ax.plot(ts, density_inverse_gaus, lw=1, c="k", label="Inverse Gaussian")
    ax.plot(ts, density_gamma, lw=1, c="k", label="Gamma")
    ax.hist(data_isi, bins=25, density=True)

    legend = ax.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0, fontsize = 9)
    legend.get_frame().set_linewidth(0.5)
    ax.set_ylabel(r"$P(I)$")
    ax.set_xlabel(r"$I$ [s]")
    plt.show()
