import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as mcolors
import scipy.interpolate as inter
import os

def print_mean_CV_to_file():
    home = os.path.expanduser("~")
    taus = []
    dCas = []
    variances = []
    rates = []
    num_spikes = []
    N = 0
    n_cha = 4
    n_clu = 10

    for i in range(50):
        for j in range(50):
            print(N)
            N += 1
            tau = np.logspace(-1, 2, 50)[j]
            dCa = np.logspace(-3, 0, 50)[i]
            #data  = np.loadtxt("/neurophysics/lukasra/Data/spikes_tau{:.2e}_j{:.2e}_Ncha4_Nclu10_Nref4_rO0.14_rC50.00_rR1.40_N0.dat".format(tau, dCa))
            data = np.loadtxt(
                home + "/CLionProjects/calcium-spikes-from-puff-phd/out/Data/spikes_adap_taua5.00e+01_e1.00e-01_tau{:.2e}_j{:.2e}_Ncha4_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(tau, dCa))
            ts, cas, jpuffs, adap = np.transpose(data)
            ISI = []
            t_tmp = 0
            for t, Ca in zip(ts, cas):
                if Ca >= 1:
                    ISI.append(t - t_tmp)
                    t_tmp = t
            if len(ISI) < 2:
                rate = 0
                variance = 0
                num_spike = 0
            else:
                rate = 1 / np.mean(ISI)
                variance = np.var(ISI)
                num_spike =  len(ISI)
            taus.append(tau)
            dCas.append(dCa)
            rates.append(rate)
            variances.append(variance)
            num_spikes.append(num_spike)

    with open(home + "/Data/Calcium/data/Ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_noadap_ncha_const.dat".format(n_clu, n_cha), "w") as outfile:
        outfile.write("# tau | dCa | 1/T | var(T) | n \n")
        for tau, dCa, r, var, n in zip(taus, dCas, rates, variances, num_spikes):
            outfile.write("{:.2e} {:.2e} {:.2f} {:.2f} {:d}\n".format(tau, dCa, r, var, n))

def plot_cv_mean_contour(text, adap):

    home = os.path.expanduser("~")
    if(adap):
        data = np.loadtxt(home + "/Data/Calcium/data/Ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_adap_ncha_const.dat".format(10, 4))
    else:
        data = np.loadtxt(
            home + "/Data/Calcium/data/Ca_mean_CV_nClu{:d}_nCha{:d}_rclose50_ncha_const.dat".format(10, 4))
    taus, dCas, rs, variances, num = np.transpose(data)

    taus = np.logspace(-1, 2, 50)
    dCas = np.logspace(-3, 0, 50)
    R = []
    R_row = []
    for r in rs:
        R_row.append(r)
        if len(R_row) == 50:
            R.append(list(R_row))
            R_row = []
    STD = []
    STD_row = []
    for var in variances:
        STD_row.append(var)
        if len(STD_row) == 50:
            STD.append(list(STD_row))
            STD_row = []
    CV = []
    CV_row = []
    for r, var, n in zip(rs, variances, num):
        if var == 0:
            CV_row.append(np.nan)
        else:
            CV_row.append(r*np.sqrt(var))
        if len(CV_row) == 50:
            CV.append(list(CV_row))
            CV_row = []
    fig, ax = plt.subplots(1,1)

    plt.scatter([taus[45], taus[45], taus[45]], [dCas[15], dCas[25], dCas[35]], c="w", zorder=5, ec="k", lw=1)
    plt.plot([taus[45], taus[45], taus[45]], [dCas[15], dCas[25], dCas[35]],c="w")
    plt.scatter([taus[30], taus[35], taus[40]], [dCas[25], dCas[25], dCas[25]],c="C7", zorder=5, ec="k", lw=1)
    plt.plot([taus[30], taus[35], taus[40]], [dCas[25], dCas[25], dCas[25]], c="C7")

    if text == "R":
        cs = ax.pcolormesh(taus, dCas, R, linewidth=0, rasterized=True, norm=mcolors.SymLogNorm(linthresh=0.01, vmin=0., vmax=1))
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
            plt.savefig(home + "/Data/Calcium/Plots/Contour_R_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Contour_R_no_adap_nClu10_nCha4.pdf", transparent=True)

    if text == "CV":
        cs = ax.pcolormesh(taus, dCas, CV, linewidth=0, rasterized=True, vmin=0., vmax=0.5)
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
            plt.savefig(home + "/Data/Calcium/Plots/Contour_CV_adap_nClu10_nCha4.pdf", transparent=True)
        else:
            plt.savefig(home + "/Data/Calcium/Plots/Contour_CV_no_adap_nClu10_nCha4.pdf", transparent=True)
    plt.show()

#print_mean_CV_to_file()
plot_cv_mean_contour("CV", adap=True)
