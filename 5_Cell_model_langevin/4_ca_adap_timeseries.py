import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    taus = np.logspace(-1, 2, 30)
    js = np.logspace(-3, 0, 30)

    tau = taus[25]
    j = js[15]
    print(tau, j)
    folder: str = "/CLionProjects/PhD/calcium_spikes_langevin/out/"
    file: str = "ca_langevin_ip1.00_taua1.00e+02_ampa2.00e-01_tau{:.2e}_j{:.2e}_N10_0.dat".format(tau, j)
    file_isi: str = "spike_times_langevin_ip1.00_taua1.00e+02_ampa2.00e-01_tau{:.2e}_j{:.2e}_N10_0.dat".format(tau, j)
    ca_data = np.loadtxt(home + folder + file)
    t, ca, jpuff, adap = np.transpose(ca_data)
    ISIs  =np.loadtxt(home + folder + file_isi)
    mean_ISI = np.mean(ISIs)
    cv_ISI = np.std(ISIs)/np.mean(ISIs)
    print(mean_ISI, cv_ISI)

    fig = plt.figure(tight_layout=True, figsize=(6, 9 / 2))
    gs = gs.GridSpec(5, 1)
    ax_adap = fig.add_subplot(gs[:2])
    ax_ca = fig.add_subplot(gs[2:])

    start = 0
    stop = 2_000
    dt = 0.1
    ax_adap.plot(t[start:stop],adap[start:stop], c="k")
    ax_ca.plot(t[start:stop], ca[start:stop], c="k")

    ax_adap.set_xlim([start*dt, stop*dt])
    ax_ca.set_xlim([start*dt, stop*dt])
    ax_ca.set_xlabel("t [s]")
    ax_ca.set_ylabel("ca")
    ax_adap.set_ylabel("a")
    ax_adap.set_xticklabels([])

    plt.show()