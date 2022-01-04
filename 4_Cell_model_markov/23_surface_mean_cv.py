import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import interp2d
from matplotlib import cm
import os
import math

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
        dR[n//50 , n%50] = r_l - r_m

    R = np.empty([50, 50])
    for n, r in enumerate(rs_m):
        R[n // 50, n % 50] = r

    CV = np.empty([50, 50])
    for n, cv in enumerate(cvs_m):
        CV[n // 50, n % 50] = cv

    f_dR = interp2d(taus, js, dR, kind='cubic')
    f_R = interp2d(taus, js, R, kind='cubic')
    taus_new = np.logspace(-1, 2, 500)
    js_new = np.logspace(-3, 0, 500)
    dR_interpolate = f_dR(taus_new, js_new)
    R_interpolate = f_R(taus_new, js_new)


    fig, ax = plt.subplots(1,1, subplot_kw={"projection": "3d"})

    if text == "R":
        taus = [math.log(x, 10) for x in taus]
        js = [math.log(x, 10) for x in js]
        TAU, JS = np.meshgrid(taus, js)

        norm = mcolors.Normalize(vmin=-0.1, vmax=0.1)
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
            plt.savefig(home + "/Data/Calcium/Plots/Surface_markov_R_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Surface_markov_R_no_adap_nClu10_nCha4.pdf", transparent=True)

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
        ax.set_ylim([-3, 0])
        ax.set_yticks([-3, -2, -1, 0])
        ax.set_yticklabels(["$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "$10^{0}$"])
        ax.set_zlabel(r"Coefficient of variation $C_V$")
        ax.set_zlim([0, 1.0])
        ax.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

        #ax.set_ylim([0.1, 10])
        #ax.set_xlim([0.01, 1])
        if (adap):
            plt.savefig(home + "/Data/Calcium/Plots/Surface_markov_CV_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Surface_markov_CV_no_adap_nClu10_nCha4.pdf", transparent=True)
    plt.show()

plot_cv_mean_contour("CV", adap=True)
