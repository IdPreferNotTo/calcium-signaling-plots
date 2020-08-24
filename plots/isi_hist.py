import numpy as np
import matplotlib.pyplot as plt
import os

def plots():
    home = os.path.expanduser("~")
    data = np.loadtxt(home + "/CLionProjects/calcium-phd/out/tauc4.0_taui10.0_alpha1.00_beta0.00_n15_m6.dat")
    isis = []
    spikes = []
    for row in data:
        if row[0] == 0:
            isi = row[1]
            isis.append(isi)
    mean = np.mean(isis)
    mu = 0.01*15*((6+1))/(2*1)
    mean_theory = 4*np.log((4*mu)/(4*mu -1.))
    var = np.var(isis)
    cv2 = var/mean**2
    cv = np.sqrt(cv2)
    plt.hist(isis, bins=50, density=True, label=r"Sim: $\langle t \rangle = {:.3f}$".format(mean) +"\n" + r"Theory: $\langle t \rangle$ = {:.3f}".format(mean_theory)  + "\n" +"$C_v =$ {:.3f}".format(cv))
    plt.legend()
    plt.show()
    return 1

if __name__ == "__main__":
    plots()