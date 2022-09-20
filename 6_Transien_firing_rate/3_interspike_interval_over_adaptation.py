import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st

if __name__ == "__main__":
    home = os.path.expanduser("~")
    folder = "/Data/calcium_spikes_markov/Data_adap"
    r0s_renewal = []
    t0s_renewal = []
    r0s_refrac = []
    t0s_refrac = []
    js = []
    taua = 100
    ampa = 0.05
    for i in range(100):
        j = 0.2*(i/100)
        js.append(j)
        r0 = np.loadtxt(home + f"/Data/calcium_spikes_theory/firing_rate/r0_ip1.00_tau2.00e+00_j{j:.2e}_N10.dat")
        r0s_renewal.append(r0)
        t0s_renewal.append(1/r0)

        #ISIs = np.loadtxt(home + folder + f"/spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}tau2.00e+00_j{j:.2e}_N10_0.dat")
        #mean_isi = np.mean(ISIs[100:])
        #r0s_refrac.append(1/mean_isi)
        #t0s_refrac.append(mean_isi)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    st.remove_top_right_axis([ax1, ax2])

    ax1.set_ylabel("$r_0$")
    ax1.plot(js, r0s_renewal, c=st.colors[0], ls=":")
    #ax1.plot(js, r0s_refrac, c=st.colors[0])
    ax1.set_xlim([0, 0.2])

    ax2.set_yscale("log")
    ax2.set_xlim([0, 0.2])
    ax2.set_xlabel(r"$P^*$")
    ax2.set_ylabel("$T_0$")
    ax2.set_ylim([1, 200])
    ax2.plot(js, t0s_renewal, c=st.colors[0], ls=":")
    #ax2.plot(js, t0s_refrac, c=st.colors[0])
    plt.savefig(home + f"/Plots/f-I_curve.pdf", transparent=True)
    plt.show()
