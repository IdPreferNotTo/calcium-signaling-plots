import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc
import default_parameters as df


if __name__ == "__main__":
    tau = 5
    p = 0.015
    M = 3
    N = 5
    K = 10
    home = os.path.expanduser("~")
    data_isi = df.load_spike_times_markov(tau, p, cer=False)
    mean = np.mean(data_isi)
    std = np.std(data_isi)
    cv = std / mean
    cv2 = cv**2
    print(data_isi.size, data_isi)

    fs = np.logspace(-3, 0, 100)
    spectrum_data = fc.power_spectrum_isis(fs, data_isi, Tmax =25 * mean)
    spectrm_inverse_gaus = fc.power_spectrum_inverse_gaussian(fs, mean, cv2)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    st.remove_top_right_axis([ax])
    ax.set_xlabel("$\omega$")
    ax.set_xscale("log")
    ax.set_ylabel("$S(\omega)$")
    #ax.set_yscale("log")
    ax.plot(fs, spectrum_data, c="C0", label="$S(\omega)$", zorder=1)
    ax.plot(fs, spectrm_inverse_gaus, lw=1, c="k", label="$S(\omega)$", zorder=0)
    ax.axhline(cv2/mean, ls=":", c="C7")
    ax.axhline(1/mean, ls=":", c="C7")

    legend = ax.legend(loc=4, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    plt.show()