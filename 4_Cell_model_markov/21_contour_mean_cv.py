import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import os

def plot_cv_mean_contour(text, adap):

    home = os.path.expanduser("~")
    home = os.path.expanduser("~")
    if(adap):
        #data = np.loadtxt(home + "/Data/Calcium/data/Ca_mean_CV_nClu{:d}_nCha{:d}_rclose50.dat".format(10, 4))
        data = np.loadtxt(home + "/Data/Calcium/data/markov_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50.dat".format(10, 4))
    else:
        #data = np.loadtxt(home + "/Data/Calcium/data/Ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_noadap.dat".format(10, 4))
        data = np.loadtxt(home + "/Data/Calcium/data/markov_ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_no_adap.dat".format(10, 4))
    taus, js, rs, cvs, num = np.transpose(data)

    taus = np.logspace(-1, 2, 50)
    js = np.logspace(-3, 0, 50)

    R = np.empty([50, 50])
    for n, r in enumerate(rs):
        R[n // 50, n % 50] = r

    CV = np.empty([50, 50])
    for n, cv in enumerate(cvs):
        if cv == 1:
            CV[n // 50, n % 50] = np.nan
        else:
            CV[n // 50, n % 50] = cv

    fig, ax = plt.subplots(1,1)

    #plt.scatter([taus[27], taus[27]], [js[30], js[47]], c="w", zorder=5, lw=1)
    #plt.plot([taus[27], taus[27]], [js[30], js[47]], c="w")

    #plt.scatter([taus[32], taus[32]], [js[25], js[47]], c="w", zorder=5, lw=1)
    #plt.plot([taus[32], taus[32]], [js[25], js[47]], c="w")

    #plt.scatter([taus[37], taus[37]], [js[20], js[47]], c="w", zorder=5, lw=1)
    #plt.plot([taus[37], taus[37]], [js[20], js[47]], c="w")
    #plt.scatter([taus[45], taus[45], taus[45]], [dCas[15], dCas[25], dCas[35]], c="w", zorder=5, ec="k", lw=1)
    #plt.plot([taus[45], taus[45], taus[45]], [dCas[15], dCas[25], dCas[35]],c="w")
    #plt.scatter([taus[30], taus[35], taus[40]], [dCas[25], dCas[25], dCas[25]],c="C7", zorder=5, ec="k", lw=1)
    #plt.plot([taus[30], taus[35], taus[40]], [dCas[25], dCas[25], dCas[25]], c="C7")

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
            plt.savefig(home + "/Data/Calcium/Plots/Contour_markov_R_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Contour_markov_R_no_adap_nClu10_nCha4.pdf", transparent=True)

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
            plt.savefig(home + "/Data/Calcium/Plots/Contour_markov_CV_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Contour_markov_CV_no_adap_nClu10_nCha4.pdf", transparent=True)
    plt.show()

plot_cv_mean_contour("CV", adap=True)
