import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    # Parameters
    tau = 1
    p = 0.060
    tau_er = 100
    eps_er = 0.1


    data_isi = df.load_spike_times_markov(tau, p, cer=True, taua=tau_er, ampa=eps_er)

    mean = np.mean(data_isi)
    std = np.std(data_isi)
    cv = std/mean
    cv2 = cv**2
    ts = np.linspace(0, 150, 501)
    p_inv_gaus = fc.inverse_gaussian_dist(ts, mean, cv2)
    p_gamma = fc.gamma_dist(ts, mean, cv2)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    st.remove_top_right_axis([ax])
    ax.set_ylabel(r"$p_{ISI}(t)$")
    ax.set_xlabel(r"$t$")

    ax.plot(ts, p_inv_gaus, lw=1, c="k", ls="--", label="Inv.\ Gaus.")
    ax.plot(ts, p_gamma, lw=1, c="k", ls=":", label="Gamma")
    ax.hist(data_isi, bins=50, color="C0", density=True)
    legend = ax.legend(fancybox=False, framealpha=1., edgecolor="k", fontsize=8)
    plt.show()