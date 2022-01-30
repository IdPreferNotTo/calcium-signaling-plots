import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import interp2d
from matplotlib import cm
import os

def plot_cv_mean_contour(text, adap):

    home = os.path.expanduser("~")
    if(adap):
        data_langevin = np.loadtxt(home + "/Data/Calcium/data/langevin_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50.dat".format(10, 4))
        data_markov = np.loadtxt(home + "/Data/Calcium/data/markov_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50.dat".format(10, 4))
    else:
        data_langevin = np.loadtxt(home + "/Data/Calcium/data/langevin_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_no_adap.dat".format(10, 4))
        data_markov = np.loadtxt(home + "/Data/Calcium/data/markov_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_no_adap.dat".format(10, 4))

    taus_l, js_l, rs_l, cvs_l, num_l = np.transpose(data_langevin)
    taus_m, js_m, rs_m, cvs_m, num_M = np.transpose(data_markov)


    taus = np.logspace(-1, 2, 50)
    js = np.logspace(-3, 0, 50)

    dR = np.empty([50, 50])
    for n, (r_l, r_m) in enumerate(zip(rs_l, rs_m)):
        dR[n//50 , n%50] = (r_l - r_m)/(r_l + r_m)

    dCV = np.empty([50, 50])
    for n , (cv_l, cv_m) in enumerate(zip(cvs_l, cvs_m)):
        dCV[n // 50, n % 50] = (cv_l -cv_m)/(cv_l + cv_m)


    fig, ax = plt.subplots(1,1)

    if text == "R":
        cs = ax.pcolormesh(taus, js, dR, linewidth=0, rasterized=True, cmap=cm.RdBu, vmin=-0.12, vmax=0.12,)#, norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))
        #ax.contour(dCas, taus, R, 5, colors="C7")
        cbar = fig.colorbar(cs)
        cbar.set_ticks([-0.1, 0, 0.1])
        cbar.set_label(r"$\Delta r =  \frac{r_{\rm langevin} - r_{\rm markov}}{r_{\rm langevin} + r_{\rm markov}}$", rotation=270, fontsize=13)
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
        cs = ax.pcolormesh(taus, js, dCV, linewidth=0, rasterized=True, vmin=-0.5, vmax=0.5, cmap=cm.RdBu)
        #ax.contour(dCas, taus, R, 5, colors="C7")
        cbar = fig.colorbar(cs)
        cbar.ax.set_yticklabels([], minor=True)
        cbar.set_label(r"$\Delta CV = \frac{CV_{\rm langevin} - CV_{\rm markov}}{CV_{\rm langevin} + CV_{\rm markov}}$", rotation=270, fontsize=13)
        ax.set_ylabel(r"$j_{\rm Blip}$")
        ax.set_xlabel(r"$\tau$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        #ax.set_ylim([0.1, 10])
        #ax.set_xlim([0.01, 1])
        if (adap):
            plt.savefig(home + "/Data/Calcium/Plots/Contour_comparison_CV_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Contour_comparison_CV_no_adap_nClu10_nCha4.pdf", transparent=True)
    plt.show()

plot_cv_mean_contour("CV", adap=True)
