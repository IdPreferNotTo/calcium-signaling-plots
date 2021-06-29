import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import interpolate
from matplotlib import cm
import os
import math

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
        CV[n // 30, n % 30] = cv


    fig, ax = plt.subplots(1,1, subplot_kw={"projection": "3d"})

    if text == "R":
        taus = [math.log(x, 10) for x in taus]
        js = [math.log(x, 10) for x in js]
        TAU, JS = np.meshgrid(taus, js)

        ax.plot_surface(TAU, JS, R, linewidth=10, cmap=cm.cividis, antialiased=True)

        ax.set_xlabel(r"Time constant $\tau$")
        ax.set_xlim([-1, 2])
        ax.set_xticks([-1, 0, 1, 2])
        ax.set_xticklabels(["$10^{-1}$", "$10^{0}$", "$10^{1}$" , "$10^{2}$"])
        ax.set_ylabel(r"Blip current $j_{\rm Blip}$")
        ax.set_ylim([-3, 0])
        ax.set_yticks([-3, -2, -1, 0])
        ax.set_yticklabels(["$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "$10^{0}$"])

        ax.set_zlabel(r"Firing rate $r$")
        #ax.set_zlim([0, 0.3])
        #ax.set_zticks([0, 0.1, 0.2, 0.3])
        #ax.set_xscale("log")
        #ax.set_yscale("log")

        if(adap):
            plt.savefig(home + "/Data/Calcium/Plots/Surface_langevin_R_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Surface_langevin_R_no_adap_nClu10_nCha4.pdf", transparent=True)

    if text == "CV":
        taus = [math.log(x, 10) for x in taus]
        js = [math.log(x, 10) for x in js]
        TAU, JS = np.meshgrid(taus, js)

        ax.plot_surface(TAU, JS, CV, linewidth=10, cmap=cm.cividis, antialiased=True)

        ax.set_xlabel(r"Time constant $\tau$")
        ax.set_xlim([-1, 2])
        ax.set_xticks([-1, 0, 1, 2])
        ax.set_xticklabels(["$10^{-1}$", "$10^{0}$", "$10^{1}$" , "$10^{2}$"])
        ax.set_ylabel(r"Blip current $j_{\rm Blip}$")
        ax.set_ylim([-3, -1])
        ax.set_yticks([-3, -2, -1])
        ax.set_yticklabels(["$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "$10^{0}$"])
        ax.set_zlabel(r"Coefficient of variation $C_V$")
        ax.set_zlim([0, 1.0])
        ax.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

        #ax.set_ylim([0.1, 10])
        #ax.set_xlim([0.01, 1])
        if (adap):
            plt.savefig(home + "/Data/Calcium/Plots/Surface_langevin_CV_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Surface_langevin_CV_no_adap_nClu10_nCha4.pdf", transparent=True)
    plt.show()

plot_cv_mean_contour("CV", adap=False)
