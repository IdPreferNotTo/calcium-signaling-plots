import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import functions as fc
import styles as st
import default_parameters as df

if __name__ == "__main__":
    Ns = np.arange(1, 100)
    Ns_mean = np.arange(1, 8)
    M = 3
    s = 1.
    Ds = []
    mus = []
    cR = 0.2
    cT = 0.5
    ci = cR
    for N in Ns:
        print(N)
        D = fc.noise_intensity_jp_single_theory(ci, N, M, s)
        mu = fc.mean_jp_single_theory(ci, N, M, s)
        Ds.append(D)
        mus.append(mu)

    home = os.path.expanduser("~")
    with open(home + f"/Data/calcium_spikes_theory/mu_D_over_N_ca{cR:.2f}.dat", "w") as f:
        f.write("#mu D N \n")
        for mu, D, N in zip(mus, Ds, Ns):
            f.write(f"{mu:.5f} {D:.5f} {N:d} \n")

    Ds_delta = []
    Ds_uniform = []
    Ds_geometric = []
    mus_delta = []
    mus_uniform = []
    mus_geometric = []
    for N_mean in Ns_mean:
        Ds_delta.append(Ds[N_mean])
        mus_delta.append(mus[N_mean])
    for N_mean in Ns_mean:
        D = 0
        mu = 0
        for n in range(1, 2*N_mean):
            p = 1/(2*N_mean-1)
            D += Ds[n]*p
            mu += mus[n]*p
        Ds_uniform.append(D)
        mus_uniform.append(mu)
    for N_mean in Ns_mean:
        D = 0
        mu = 0
        for n in range(1, 10*N_mean):
            p = np.power((1 - 1/N_mean), n-1)/N_mean
            D += Ds[n]*p
            mu += mus[n]*p
        Ds_geometric.append(D)
        mus_geometric.append(mu)


    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(9, 3))
    gs = gridspec.GridSpec(1, 3)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    st.remove_top_right_axis([ax1, ax2, ax3])
    color = st.Colors

    ax1.set_ylabel(r"$\mu_{\rm puff}/\langle N \rangle$")
    ax1.set_xlabel(r"$\langle N \rangle$")
    ax1.set_xticks(Ns_mean)
    ax1.plot(Ns_mean, [mu/n for mu, n in zip(mus_delta, Ns_mean)], color=color.palette[1], label="delta")
    ax1.plot(Ns_mean, [mu/n for mu, n in zip(mus_uniform, Ns_mean)], color=color.palette[3], label="uniform")
    ax1.plot(Ns_mean, [mu/n for mu, n in zip(mus_geometric, Ns_mean)], color=color.palette[4], label="geometric")
    ax1.legend(fancybox=False, framealpha=1., fontsize=11)
    ax1.set_ylim([0, 0.07])

    ax2.set_ylabel(r"$D_{\rm puff}/\langle N \rangle$")
    ax2.set_xlabel(r"$\langle N \rangle$")
    ax2.set_xticks(Ns_mean)
    ax2.plot(Ns_mean, [D/n for D, n in zip(Ds_delta, Ns_mean)], color=color.palette[1])
    ax2.plot(Ns_mean, [D/n for D, n in zip(Ds_uniform, Ns_mean)], color=color.palette[3])
    ax2.plot(Ns_mean, [D/n for D, n in zip(Ds_geometric, Ns_mean)], color=color.palette[4])
    ax2.set_ylim([0, 0.07])

    ax3.set_ylabel(r"$D_{\rm puff} / \mu_{\rm puff}$")
    ax3.set_xlabel(r"$\langle N \rangle$")
    ax3.set_xticks(Ns_mean)
    ax3.plot(Ns_mean, [D / mu for D, mu, n in zip(Ds_delta, mus_delta, Ns_mean)], color=color.palette[1])
    ax3.plot(Ns_mean, [D / mu for D, mu, n in zip(Ds_uniform, mus_uniform, Ns_mean)], color=color.palette[3])
    ax3.plot(Ns_mean, [D / mu for D, mu, n in zip(Ds_geometric, mus_geometric, Ns_mean)], color=color.palette[4])
    ax3.set_ylim([0, 1.0])
    plt.show()