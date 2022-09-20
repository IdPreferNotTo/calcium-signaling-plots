import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import os

import styles as st
import functions as fc

def transient_func2(i, T0, T8, tau):
    return T0*np.exp(-i/tau) + T8*(1 - np.exp(-i/tau))

def get_r0_theory(tau, current):
    cas_theory = np.linspace(0.30, 1, 1401)
    dca = cas_theory[1] - cas_theory[0]
    p0s_theo_ca = []
    Ds_theo_ca = []
    integral = 0
    n_clu = 10
    for ca in reversed(cas_theory[1:]):
        print(f"{ca:.3f}")
        h = fc.h_func(ca, tau, current, n_clu, 5, 4, 1)
        d = fc.d_func(ca, current, n_clu, 5, 4, 1)
        if ca == 1:
            integral += 0
        elif ca >= 0.33:
            integral += np.exp(-h) * dca
        Ds_theo_ca.append(d)
        p0s_theo_ca.append(integral * np.exp(h) / d)
    norm = np.sum(p0s_theo_ca) * dca
    return 1 / norm

class GradientDesecent():
    def __init__(self, function):
        self.dj = 0.01
        self.dtau = 0.1
        self.function = function

    def opt(self, tau, j, r0_target):
        r0 = self.function(tau, current)

        grad_tau_plus = (r0 - self.function(tau + self.dtau, j))/self.dtau
        grad_j_plus = (r0 - self.function(tau, j + self.dj))/self.dj
        grad_tau_minus = (r0 - self.function(tau + self.dtau, j)) / self.dtau
        grad_j_minus = (r0 - self.function(tau, j + self.dj)) / self.dj
        grad_tau = (grad_tau_plus + grad_tau_minus)/2
        grad_j = (grad_j_plus + grad_j_minus)/2



if __name__ == "__main__":
    home = os.path.expanduser("~")
    file_str = home + "/Data/calcium_spikes_experimental/Spikes/HEK/HEK2_bapta_ratio.dat"
    data = np.loadtxt(file_str)
    n = len(data[0])
    j = 12
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    print(idxs)
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx - i), axis=0)
    caT = 0.4
    stat_cas = []
    times = [x[0] for x in data]
    cas = [x[j] for x in data]
    spiking = False
    spike_times = []
    spike_end_times = []
    for idx, (t, y) in enumerate(zip(times, cas)):
        if t < 3500:
            if y > caT and not spiking:
                spike_times.append(t)
                spiking = True
            if y < caT and spiking:
                spiking = False

    ISIs = []
    for t1, t2 in zip(spike_times[:-1], spike_times[1:]):
        ISIs.append(t2 - t1)

    nr_ISIs = len(ISIs)
    index_ISIs = np.arange(nr_ISIs)
    popt, pcov = curve_fit(transient_func2, index_ISIs, ISIs, p0=(100, 150, 2))

    tau = 5.
    current = 0.036
    r0theo = get_r0_theory(tau, current)
    r0data = 1/popt[0]
    print(r0data, r0theo)

