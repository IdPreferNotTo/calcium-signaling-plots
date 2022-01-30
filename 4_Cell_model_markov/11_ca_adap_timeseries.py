import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
import os
from functions import *

import styles as st

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(8, 4))
    gs = gridspec.GridSpec(4, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1:3, 0])
    ax3 = fig.add_subplot(gs[3, 0])
    ax4 = fig.add_subplot(gs[0:2, 1])
    st.remove_top_right_axis([ax1, ax2, ax3, ax4])

    tau = 0.2
    j = 0.501
    taua = 500
    ampa = 0.05
    home = os.path.expanduser("~")
    print("Load Data...")
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/Data_adap_fix_adap_para/"
    file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    file_spike =  f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
    data = np.loadtxt(home + folder + file)
    isis  = np.loadtxt(home + folder + file_spike)
    mean_isi = np.mean(isis)
    cv_isi = np.std(isis)/mean_isi
    print(mean_isi, cv_isi)
    ts, cas, jpuffs, adaps = np.transpose(data)
    print("done")

    spike_times = []
    cR = 0.33
    cT = 1.0
    t_isi = []
    ca_isi = []
    jpuff_isi = []
    adap_isi= []
    for t, ca, jpuff, adap, adap_after in zip(ts[:-1], cas[:-1], jpuffs[:-1], adaps[:-1], adaps[1:]):
        if t > 500:
            break
        t_isi.append(t)
        ca_isi.append(ca)
        jpuff_isi.append(jpuff)
        adap_isi.append(adap)
        if ca == 1:
            spike_times.append(t)
            ax1.plot(t_isi, adap_isi, c="C0")
            ax1.plot([t, t], [adap, adap_after], ls=":", c="C0")

            ax2.arrow(x=t, y=1, dx=0, dy=1, color="C0", length_includes_head = True, head_width = 7.5, head_length=0.1)
            ax2.plot(t_isi, ca_isi, c="C0")
            ax2.plot([t, t], [cR, cT], ls=":", c="C0")

            ax3.plot(t_isi, jpuff_isi, c="C0")
            t_isi.clear()
            ca_isi.clear()
            jpuff_isi.clear()
            adap_isi.clear()


    ax1.set_ylabel(r"$c_{ER}$")
    ax1.set_xlim([0, 500])
    ax1.set_ylim([0.5, 1])

    ax2.set_ylabel(r"$c_i$")
    ax2.set_ylim([0.8*cR, 2.1*cT])
    ax2.set_yticks([cR, cT])
    ax2.set_yticklabels(["$c_R$", "$c_T$"])
    ax2.set_xlim([0, 500])

    ax3.set_xlabel("$t$ / s")
    ax3.set_xlim([0, 500])
    ax3.set_ylabel(r"$j_{\rm puff}$")

    ax4.hist(isis, density=True)
    mean_isi = np.mean(isis)
    cv_isi = np.std(isis)/mean_isi
    xs_exp = np.linspace(0, 3*mean_isi)
    ys_exp = 1/(mean_isi)*np.exp(-xs_exp/mean_isi)
    xs_ing, ys_ing = inverse_gaussian(mean_isi, cv_isi)
    ax4.plot(xs_exp, ys_exp)
    ax4.plot(xs_ing, ys_ing)
    plt.savefig(home + f"/Data/Calcium/Plots/7_markov_ca_adap_timeseries_tau{tau:.2e}j{j:.2e}.pdf", transparent=True)
    plt.show()
