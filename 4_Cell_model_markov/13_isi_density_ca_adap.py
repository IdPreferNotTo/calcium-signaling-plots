import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from functions import *

if __name__ == "__main__":
    home = os.path.expanduser("~")
    # Parameters
    tau = 2.81
    j = 0.0728
    taua = 100
    ampa = 0.2
    home = os.path.expanduser("~")

    folder = "/Data/calcium_spikes_markov/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    isis = np.loadtxt(home + folder + file_spike)

    mean_isi = np.mean(isis)
    std_isi = np.std(isis)
    cv_isi = std_isi/mean_isi
    cv2_isi = cv_isi**2
    ts = np.linspace(0, 150, 501)
    inv_gaus = []
    for t in ts:
        p = np.sqrt(mean_isi/(2*np.pi*cv2_isi*(t**3)))*np.exp(-(t - mean_isi)**2/(2*mean_isi*cv2_isi*t))
        inv_gaus.append(p)

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    remove_top_right_axis([ax])
    ax.set_ylabel(r"$P(I)$")
    ax.set_xlabel(r"$I$ [s]")

    ax.plot(ts, inv_gaus, lw=1, c="k", label="Inverse Gaussian")
    ax.hist(isis, bins=50, color="C0", density=True)
    plt.show()