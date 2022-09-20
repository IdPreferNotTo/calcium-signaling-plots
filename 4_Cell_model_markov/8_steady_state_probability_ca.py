import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import styles as st
import functions as fc


if __name__ == "__main__":
    home = os.path.expanduser("~")
    # Parameters
    taus = np.logspace(-1, 2, 50)
    jcas = np.logspace(-3, 0, 50)
    tau = taus[30]
    jca = jcas[25]
    m = 4
    n = 5
    N = 10
    print(tau, jca)

    # Load Data
    folder = "/Data/calcium_spikes_markov/Data_no_adap/"
    file_ca = "ca_markov_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, jca, N)
    file_spikes = "spike_times_markov_ip1.00_tau{:.2e}_j{:.2e}_N{:d}_0.dat".format(tau, jca, N)
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
        f_theo_ca.append(fc.f_func(ca, tau, jca, N, n, m, 1))
        sqrt2D_theo_ca.append(np.sqrt(2 * fc.d_func(ca, jca, N, n, m, 1)))


    cas_theory = np.linspace(0.30, 1, 1401)
    dca = cas_theory[1] - cas_theory[0]
    p0s_theo_ca = []
    Ds_theo_ca = []
    I2s_theo_ca = []
    integral = 0

    for ca in reversed(cas_theory[1:]):
        print(f"{ca:.3f}")
        h = fc.h_func(ca, tau, jca, N, n, m, 1)
        d = fc.d_func(ca, jca, N, n, m, 1)
        if ca == 1:
            integral += 0
        elif ca >= 0.33:
            integral += np.exp(-h)*dca
        Ds_theo_ca.append(d)
        p0s_theo_ca.append(integral * np.exp(h) / d)

    integral = 0
    for ca in cas_theory[1:]:
        print(f"{ca:.3f}")
        h = fc.h_func(ca, tau, jca, N, n, m, 1)
        d = fc.d_func(ca, jca, N, n, m, 1)
        integral += (np.exp(h) / d)*dca
        I2s_theo_ca.append(np.power(integral*np.exp(-h),2))

    norm = np.sum(p0s_theo_ca) * dca
    r0theo = 1/norm
    ISI_mean_theo  =1/r0theo
    ISI_mean = np.mean(isi_data)
    cvsim = np.std(isi_data)/np.mean(isi_data)
    p0s_theo_ca = [x / norm for x in p0s_theo_ca]

    cv2_theo = 0
    for I2s, D, P0 in zip(I2s_theo_ca, Ds_theo_ca, p0s_theo_ca):
        cv2_theo += 2*r0theo*I2s*D*P0*dca
    cv_theo = np.sqrt(cv2_theo)
    print(cv_theo)

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    st.remove_top_right_axis([ax])


    label_ca_sim = "Simulation \n" + rf"$\langle T \rangle$ = {ISI_mean:.1e}" + "\n" + rf"$C_{{V}}$ = {cvsim:.1e}"
    label_ca_theo = "Theory \n" + rf"$\langle T \rangle$ = {ISI_mean_theo:.1e}" + "\n" + rf"$C_{{V}}$ = {cv_theo:.1e}"
    ax.plot(cas_theory[1:], p0s_theo_ca[::-1], lw=1, c="k", label=label_ca_theo)
    ax.hist(cas, bins=50, density=True, alpha=.6, color="C0",
            label=label_ca_sim)
    legend = ax.legend(loc=1, fancybox=False, edgecolor="k", framealpha=1.0)
    legend.get_frame().set_linewidth(0.5)
    ax.set_ylabel(r"$P_0([\rm{Ca}^{2+}]_i)$")
    ax.set_xlabel(r"$[\rm{Ca}^{2+}]_i$ [a.u.]")
    plt.savefig(home + "/Data/Calcium/Plots/ca_steady_state_probability.pdf", transparent=True)
    plt.show()
