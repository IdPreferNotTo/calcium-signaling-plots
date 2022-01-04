import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    folder = home + "/Data/Calcium/data"
    file = "/markov_adap_parameters_tau1.05e+01_j1.46e-02.dat"
    data = np.loadtxt(folder + file) #tau, j, T, CV, n

    d = 100
    T0 = 76
    tau_eff0 = 8.2
    for set in data:
        T = set[3]
        tau_eff = set[4]
        dT = abs(T-T0)/T0
        dtau = abs(tau_eff-tau_eff0)/tau_eff0
        if (dT + dtau)/2. < d:
            taua0 = set[0]
            ampa0 = set[1]
            d = (dT + dtau)/2.

    print(taua0, ampa0)
