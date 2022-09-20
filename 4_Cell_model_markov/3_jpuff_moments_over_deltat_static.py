import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


import functions as fc
import styles as st

if __name__ == "__main__":
    ca_fix = 0.33
    ca_res = 0.33
    N = 10
    n = 5
    jca = 1
    tau = 1
    home = os.path.expanduser("~")
    folder = "/Data/calcium_spikes_markov/ca_fix/"
    file = "ca_markov_cafix{:.2f}_ip1.00_tau1.00e+00_j1.00e+00_N10_0.dat".format(ca_fix)
    data = np.loadtxt(home + folder + file)
    ts, cas, jpuffs, adaps = np.transpose(data)

    ts_clear = []
    cas_clear = []
    jpuffs_clear = []
    mean_jpuffs_clear = []
    for t, ca, jpuff in zip(ts, cas, jpuffs):
        if ca != 1:
            ts_clear.append(t)
            cas_clear.append(ca)
            jpuffs_clear.append(jpuff)
    mean = np.mean(jpuffs_clear)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax0 = fig.add_subplot(gs[0, 0])
    st.remove_top_right_axis([ax0])

    dt0 = 0.001
    dts = []
    mu2s = []
    mu3s = []
    mu4s = []

    factors = np.arange(1, 100)
    for f in factors:
        cg_jpuffs = fc.coarse_grain_list(jpuffs_clear, f)
        mu2 = fc.moments(cg_jpuffs, 2)
        mu2s.append(mu2)
        mu3 = fc.moments(cg_jpuffs, 3)
        mu3s.append(mu3)
        mu4 = fc.moments(cg_jpuffs, 4)
        mu4s.append(mu4)
        dts.append(f * dt0)

    ax0.plot(dts, mu2s, label=r"$\gamma_2 = \langle (x - \mu)^2  \rangle$")
    ax0.plot(dts, mu3s, label=r"$\gamma_3 = \langle (x - \mu)^3  \rangle$")
    ax0.plot(dts, mu4s, label=r"$\gamma_4 = \langle (x - \mu)^4  \rangle$")
    legend = ax0.legend(loc=3, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    # ax0.set_ylim([0.001, 1])
    ax0.set_ylabel("$\mu_x$")
    ax0.set_yscale("log")
    ax0.set_xlabel("$\Delta t$")
    ax0.set_xscale("log")
    #plt.savefig(home + "/Data/Calcium/Plots/puff_current.pdf".format(ca_fix, tau, jca), transparent=True)
    plt.show()
