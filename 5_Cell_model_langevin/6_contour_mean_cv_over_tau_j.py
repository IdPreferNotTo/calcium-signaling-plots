import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import os

def plot_cv_mean_contour(text, adap):

    home = os.path.expanduser("~")
    if(adap):
        #data = np.loadtxt(home + "/Data/Calcium/data/Ca_mean_CV_nClu{:d}_nCha{:d}_rclose50.dat".format(10, 4))
        data = np.loadtxt(home + "/Data/Calcium/data/langevin_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50.dat".format(10, 4))
    else:
        #data = np.loadtxt(home + "/Data/Calcium/data/Ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_noadap.dat".format(10, 4))
        data = np.loadtxt(home + "/Data/Calcium/data/langevin_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_no_adap.dat".format(10, 4))
    taus, js, rs, cvs, num = np.transpose(data)

    taus = np.logspace(-1, 2, 30)
    js = np.logspace(-3, 0, 30)

    R = np.empty([30, 30])
    for n, r in enumerate(rs):
        R[n // 30, n % 30] = r

    CV = np.empty([30, 30])
    for n, cv in enumerate(cvs):
        if cv == 1:
            CV[n // 30, n % 30] = np.nan
        else:
            CV[n // 30, n % 30] = cv

    fig, ax = plt.subplots(1,1)

    if text == "R":
        cs = ax.pcolormesh(taus, js, R, linewidth=0, rasterized=True, cmap=cm.cividis)#, norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))
        #ax.contour(dCas, taus, R, 5, colors="C7")
        cbar = fig.colorbar(cs)
        cbar.set_label('$r=1/T$', rotation=270)
        ax.set_ylabel(r"$j_{\rm Blip}$")
        ax.set_xlabel(r"$\tau$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        #ax.set_ylim([0.1, 10])
        #ax.set_xlim([0.01, 1])
        if(adap):
            plt.savefig(home + "/Data/Calcium/Plots/Contour_langevin_R_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Contour_langevin_R_no_adap_nClu10_nCha4.pdf", transparent=True)

    if text == "CV":
        cs = ax.pcolormesh(taus, js, CV, linewidth=0, rasterized=True, vmin=0., vmax=0.5, cmap=cm.cividis)
        #ax.contour(dCas, taus, R, 5, colors="C7")
        cbar = fig.colorbar(cs)
        cbar.ax.set_yticklabels([], minor=True)
        cbar.set_label('$CV$', rotation=270)
        ax.set_ylabel(r"$j_{\rm Blip}$")
        ax.set_xlabel(r"$\tau$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        #ax.set_ylim([0.1, 10])
        #ax.set_xlim([0.01, 1])
        if (adap):
            plt.savefig(home + "/Data/Calcium/Plots/Contour_langevin_CV_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Contour_langevin_CV_no_adap_nClu10_nCha4.pdf", transparent=True)
    plt.show()

plot_cv_mean_contour("R", adap=True)
