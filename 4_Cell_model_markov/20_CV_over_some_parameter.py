import os
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs


if __name__ == "__main__":
    home = os.path.expanduser("~")
    taus = np.logspace(-1, 2, 30)
    js = np.logspace(-3, 0, 30)
    means = []
    stds = []
    means_no_adap = []
    stds_no_adap = []
    for i in range(7):
        tau = taus[17+i]
        j = js[17+i]
        taua = 100
        ampa = 0.2
        folder = "/Data/calcium_spikes_markov/"
        file = f"ca_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
        file_spike = f"spike_times_markov_ip1.00_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{j:.2e}_N10_0.dat"
        data = np.loadtxt(home + folder + file)
        isis = np.loadtxt(home + folder + file_spike)
        mean_isi = np.mean(isis)
        std_isi = np.std(isis)
        cv_isi = std_isi / mean_isi
        print(mean_isi, cv_isi)
        means.append(mean_isi)
        stds.append(std_isi)


    fig = plt.figure(tight_layout=True, figsize=(6*3/4, 9*3/(4*2)))
    gs = gs.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:])
    ax.scatter(means, stds)

    p0, p1 = poly.polyfit(means, stds, 1)
    xs = np.linspace(0, max(means))
    ys = [p0 + p1*x for x in xs]
    ax.plot(xs, ys, c="C7", label = "slope = {:.2f}".format(p1))

    #x.scatter(means_no_adap, stds_no_adap)
    #p0, p1 = poly.polyfit(means_no_adap, stds_no_adap, 1)
    #xs = np.linspace(0, max(means))
    #ys = [p0 + p1*x for x in xs]
    #ax.plot(xs, ys, c="C1", label="slope = {:.2f}".format(p1))
    ax.legend()
    ax.set_xlabel(r"$\mu = \langle {\rm ISI} \rangle$")
    ax.set_ylabel(r"$\sigma = \langle (\Delta {\rm ISI})^2 \rangle$")
    #plt.savefig(home + "/Data/Calcium/Plots/CV_jblip_tau_const_ratio_1.pdf", transparent=True)
    plt.show()

