import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import os

def plot_cv_mean_contour(text, adap):

    home = os.path.expanduser("~")
    if(adap):
        data = np.loadtxt(home + "/Data/Calcium/data/langevin_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50.dat".format(10, 4))
    else:
        data = np.loadtxt(home + "/Data/calcium_spikes_theory/langevin_strat_ca_mean_CV_K{:d}_N{:d}_no_adap.dat".format(10, 5))

    taus, js, isis, cvs, num = np.transpose(data)

    N = 61
    taus = np.logspace(-1, 2, N)
    ps = np.logspace(-3, 0, N)


    ISIs = np.empty([N, N])
    for n, I in enumerate(isis):
        ISIs[n // N, n % N] = I

    CV = np.empty([N, N])
    for n, cv in enumerate(cvs):
        if cv == 1:
            CV[n // N, n % N] = np.nan
        else:
            CV[n // N, n % N] = cv

    fig, ax = plt.subplots(1,1)

    if text == "R":
        cs = ax.pcolormesh(taus, ps, ISIs, linewidth=0, rasterized=True, cmap=cm.cividis, norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=1., vmax=1000))
        cbar = fig.colorbar(cs)
        cbar.set_label(r'$\langle I \rangle$', rotation=270)
        ax.set_ylabel(r"$j_{\rm Blip}$")
        ax.set_xlabel(r"$\tau$")
        ax.set_xscale("log")
        ax.set_yscale("log")

    if text == "CV":
        cs = ax.pcolormesh(taus, ps, CV, linewidth=0, rasterized=True, cmap=cm.cividis, norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=0.1, vmax=1))
        cbar = fig.colorbar(cs)
        cbar.ax.set_yticklabels([], minor=True)
        cbar.set_label('$CV$', rotation=270)
        ax.set_ylabel(r"$j_{\rm Blip}$")
        ax.set_xlabel(r"$\tau$")
        ax.set_xscale("log")
        ax.set_yscale("log")
    plt.show()

plot_cv_mean_contour("R", adap=False)
