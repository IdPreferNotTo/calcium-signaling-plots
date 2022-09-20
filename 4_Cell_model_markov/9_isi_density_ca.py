import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

import functions as fc
import styles as st





if __name__ == "__main__":
    home = os.path.expanduser("~")
    # Parameters
    taus = np.logspace(-1, 2, 50)
    jcas = np.logspace(-3, 0, 50)
    tau = 1.68
    jca = 0.0596
    m = 4
    n = 5
    N = 10
    print(tau, jca)
    # Load Data
    folder = "/Data/calcium_spikes_markov/Data_no_adap/"
    file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{jca:.2e}_N{N:d}_0.dat"
    isi_data = np.loadtxt(home + folder + file_spikes)

    mean_isi = np.mean(isi_data)
    std_isi = np.std(isi_data)
    cv_isi = std_isi/mean_isi
    cv2_isi = cv_isi**2
    ts = np.linspace(0, 50, 501)
    inv_gaus = []
    for t in ts:
        p = np.sqrt(mean_isi/(2*np.pi*cv2_isi*(t**3)))*np.exp(-(t - mean_isi)**2/(2*mean_isi*cv2_isi*t))
        inv_gaus.append(p)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    st.remove_top_right_axis([ax])


    ax.plot(ts, inv_gaus, lw=1, c="k", label="Inverse Gaussian")
    ax.hist(isi_data, bins=50, density=True, alpha=.6, color="C0", label=f"$P(I)$ \n $<T> =$ {mean_isi:.1f} \n $CV_T = {cv_isi:.2f}$")

    legend = ax.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    ax.set_ylabel(r"$P(I)$")
    ax.set_xlabel(r"$I$ [s]")
    #plt.savefig(home + "/Data/Calcium/Plots/ca_isi_probability.pdf", transparent=True)
    plt.show()
