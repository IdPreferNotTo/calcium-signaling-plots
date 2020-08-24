import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import os

def plots():
    home = os.path.expanduser("~")
    data = np.loadtxt(home + "/CLionProjects/calcium-phd/out/tauc4.0_taui10.0_alpha0.00_beta0.00_n15_m6.dat")
    ipis_by_site = []
    amplitude_by_site = []
    for i in range(15):
        ipis_by_site.append([])
        amplitude_by_site.append([])
    for row in data:
        if row[0] != 0:
            index = int(row[0])-1
            ipi_by_site = row[1]
            puff_amplitude = row[2]
            amplitude_by_site[index].append(puff_amplitude)
            ipis_by_site[index].append(ipi_by_site)
        else:
            continue

    site = 1
    mean_ipi = np.mean(ipis_by_site[site])
    channels_per_site =[]
    for a in amplitude_by_site:
        channels_per_site.append(np.amax(a)/0.01)
    N = np.mean(channels_per_site)

    fig = plt.figure()
    for c, site in zip( ["C0"], [1]):
        mean_ipi_theory = 1*(6/1)*(1/channels_per_site[site])
        var_ipi_theory = 2*0.5*((6/1)*(1/channels_per_site[site]))**3
        Cv2 = var_ipi_theory/(mean_ipi_theory**2)
        ts = np.linspace(0.01, 4*mean_ipi, 100)
        igaussian = []
        for t in ts:
            value = np.sqrt(mean_ipi_theory/(2*np.pi*Cv2*(t**3)))*np.exp(-(((t-mean_ipi_theory)**2)/(2*Cv2*mean_ipi_theory*t)))
            igaussian.append(value)
        plt.hist(ipis_by_site[site], bins=100, density=True, color=c, alpha=0.5)
        plt.plot(ts, igaussian, color=c)
    plt.xlim(0, 10)
    plt.show()
    return 1

if __name__ == "__main__":
    plots()