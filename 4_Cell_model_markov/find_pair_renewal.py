import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    folder = home + "/Data/Calcium/data"
    file = "/markov_ca_mean_CV_N10_n5_no_adap.dat"
    data = np.loadtxt(folder + file) #tau, j, T, CV, n

    d = 100
    tau0 = 0
    j0 = 0
    T0 = 35
    CV0 = 0.1
    for set in data:
        T = set[2]
        CV = set[3]
        dT = abs(T-T0)/T0
        dCV = abs(CV-CV0)/CV0
        if (dT + dCV)/2. < d:
            tau0 = set[0]
            j0 = set[1]
            d = (dT + dCV)/2.

    print(tau0, j0)
