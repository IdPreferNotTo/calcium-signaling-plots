import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from functions import *


if __name__ == "__main__":
    home = os.path.expanduser("~")
    # Parameters
    tau = 2.81
    j = 0.0728 #0.0574
    n_cha = 4
    n_clu = 10
    print(tau, j)
    # Load Data
    folder = "/CLionProjects/PhD/calcium_spikes_markov/out/"
    file_ca = "ca_markov_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu)
    file_spikes = "spike_times_markov_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, j, n_clu)
    ca_data = np.loadtxt(home + folder + file_ca)
    ca_data = [set for set in ca_data if set[1] != 1.]
    isi_data = np.loadtxt(home + folder + file_spikes)
    ts, cas, jpuffs, adaps = np.transpose(ca_data)


    # f is the whole deterministic dynamics including leak and mean puff current
    cas_plot = np.linspace(0, 1, 101, endpoint=True)

    p0s_sim_ca = [[] for _ in range(100)]
    fs_sim_ca =  []
    sqrt2ds_sim_ca = []

    for ca, jpuff in zip(cas, jpuffs):
        for i, (ca_min, ca_max) in enumerate(zip(cas_plot[:-1], cas_plot[1:])):
            if ca_min < ca and ca < ca_max:
                p0s_sim_ca[i].append(jpuff)

    for ca, p0_sim_ca in zip(cas_plot, p0s_sim_ca):
        if len(p0_sim_ca) == 0:
            fs_sim_ca.append(np.nan)
            sqrt2ds_sim_ca.append(np.nan)
        else:
            fs_sim_ca.append(np.mean(p0_sim_ca) - (ca - 0.33) / tau)
            sqrt2ds_sim_ca.append(np.std(p0_sim_ca))

    # Calculate f(x) and D(x) from theory. Convention dx/dt = f(x) + \sqrt{2D(x)}xi(t)
    # Hence f(x) = j_leak + <j_puff> and D(x) = j^2 * N * D_puff(x)

    sqrt2D_theo_ca = []
    f_theo_ca = []
    for ca in cas_plot:
        f_theo_ca.append(f_func(ca, tau, j, n_clu))
        sqrt2D_theo_ca.append(np.sqrt(2*d_func(ca, j, n_clu)))


    cas_theory = np.linspace(0.30, 1, 1001)
    dca = cas_theory[1] - cas_theory[0]
    p0s_theo_ca = []
    Ds_theo_ca = []
    I2s_theo_ca = []
    integral = 0

    for ca in reversed(cas_theory[1:]):
        print(ca)
        h = h_func(ca, tau, j, n_clu)
        d = d_func(ca, j, n_clu)
        if ca == 1:
            integral += 0
        elif ca >= 0.33:
            integral += np.exp(-h)*dca
        Ds_theo_ca.append(d)
        p0s_theo_ca.append(integral * np.exp(h) / d)

    integral = 0
    for ca in cas_theory[1:]:
        print(ca)
        h = h_func(ca, tau, j, n_clu)
        d = d_func(ca, j, n_clu)
        integral += (np.exp(h) / d)*dca
        I2s_theo_ca.append(np.power(integral*np.exp(-h),2))

    #mean_T = 0
    #integral = 0
    #for ca in cas_theory[1:]:
    #    print(ca)
    #    h = h_func(ca, tau, j, n_clu)
    #    d = d_func(ca, j, n_clu)
    #    if ca < 0.33:
    #        integral += (np.exp(h)/d)*dca
    #    else:
    #        mean_T += np.exp(-h)*integral*dca
    #        integral += (np.exp(h)/d)*dca

    #print(1/mean_T)

    norm = np.sum(p0s_theo_ca) * dca
    r0theo = 1/norm
    r0sim = 1/np.mean(isi_data)
    cvsim = np.std(isi_data)/np.mean(isi_data)
    print(r0theo)
    print(r0sim)
    p0s_theo_ca = [x / norm for x in p0s_theo_ca]

    print(I2s_theo_ca)
    print(Ds_theo_ca)
    print(p0s_theo_ca)
    cv2_theo = 0
    for I2s, D, P0 in zip(I2s_theo_ca, Ds_theo_ca, p0s_theo_ca):
        cv2_theo += 2*r0theo*I2s*D*P0*dca
    cv_theo = np.sqrt(cv2_theo)
    print(cv_theo)

    set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    #axins = inset_axes(ax, width="40%", height="40%", loc=1)
    remove_top_right_axis([ax])

    #axins.plot(cas_plot[:-1], fs_sim_ca, c="k", label="f(x)")
    #axins.plot(cas_plot, f_theo_ca, c="k", alpha=0.7, ls="--")
    #axins.plot(cas_plot[:-1], sqrt2ds_sim_ca, label="$\sqrt{2D(x)}$", c="C0", )
    #axins.plot(cas_plot, sqrt2D_theo_ca, c="C0", alpha=0.7, ls="--")
    #axins.legend()

    label_ca_sim = "Simulation \n" + rf"$r_{{0}}$ = {r0sim:.1e}" + "\n" + rf"$C_{{V}}$ = {cvsim:.1e}"
    label_ca_theo = "Theory \n" + rf"$r_{{0}}$ = {r0theo:.1e}" + "\n" + rf"$C_{{V}}$ = {cv_theo:.1e}"
    ax.plot(cas_theory[1:], p0s_theo_ca[::-1], lw=1, c="k", label=label_ca_theo)
    ax.hist(cas, bins=50, density=True, alpha=.6, color="C0",
            label=label_ca_sim)
    legend = ax.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    ax.set_ylabel(r"$P_0([\rm{Ca}^{2+}]_i)$")
    ax.set_xlabel(r"$[\rm{Ca}^{2+}]_i$ [a.u.]")
    plt.savefig(home + "/Data/Calcium/Plots/7_markov_ca_steady_state_probability.pdf", transparent=True)
    plt.show()
