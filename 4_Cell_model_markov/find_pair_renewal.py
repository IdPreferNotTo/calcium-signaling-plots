import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    folder = home + "/Data/calcium_spikes_theory"
    file = "/markov_ca_mean_CV_K10_N5_no_adap.dat"
    data = np.loadtxt(folder + file) #tau, j, T, CV, n

    d = 100
    tau0 = 0
    j0 = 0
    T0 = 35.
    CV0 = 0.10

    pairs = []
    parameters_pair = []
    for set in data:
        T = set[2]
        CV = set[3]
        if 0.18 < CV < 0.22:
            pairs.append([set[2], set[3]])
            parameters_pair.append([set[0], set[1]])
    print(pairs)
    print(parameters_pair)
    
    for set in data:
        T = set[2]
        CV = set[3]
        dT = abs(T-T0)/T0
        dCV = abs(CV-CV0)/CV0
        if (dT + dCV)/2. < d:
            tau0 = set[0]
            j0 = set[1]
            T_found = T
            CV_found = CV
            d = (dT + dCV)/2.

    print(T_found, CV_found)
    print("tau:", tau0, "p:", j0)
