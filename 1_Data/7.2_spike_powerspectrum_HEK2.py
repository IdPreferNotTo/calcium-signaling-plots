import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from functions import *


home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
data = np.loadtxt(file_str)

ISIs = []
MEANs = []
STDs = []
CVs = []
n = len(data[0])
correlations = []
for j in range(1, n):
    ts = [x[0] for x in data]
    cas = [x[j] for x in data]
    spiking: bool = False
    t_tmp: float = 0
    spike_times = []
    if j == 2:
        for t, ca in zip(ts, cas):
            if t > 500 and t < 3700:
                if ca > 0.35 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.35 and spiking:
                    spiking = False
    else:
        for t, ca in zip(ts, cas):
            if t > 500 and t < 3700:
                if ca > 0.4 and not spiking:
                    spike_times.append(t)
                    spiking = True
                if ca < 0.4 and spiking:
                    spiking = False

    isis = [t2 - t1 for t2, t1 in zip(spike_times[1:], spike_times[:-1])]
    isi_var = k_corr(isis, isis, 0)
    isi_covar = k_corr(isis, isis, 1)
    correlations.append(isi_covar/isi_var)

    ISIs.append([t for t in isis])
    n = len(isis)
    mean = np.mean(isis)
    MEANs.append(mean)
    std = np.std(isis)
    STDs.append(std)
    cv = std / mean
    cv2 = np.power(cv, 2)
    CVs.append(cv)
coef = np.polyfit(MEANs, STDs, 1)
t_ref = -coef[1]/coef[0]

set_default_plot_style()
fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])


remove_top_right_axis([ax])
ISIs_no_ref = []
for isis in ISIs:
    isis_no_ref = []
    for isi in isis:
        isis_no_ref.append(isi - t_ref)
    ISIs_no_ref.append([t/np.mean(isis_no_ref) for t in isis_no_ref])
ISIs_flat = [item for sublist in ISIs_no_ref for item in sublist]

ISIs_chunks = []

chunks = []
t = 0
Tmax = 25
for isi in ISIs_flat:
    t += isi
    chunks.append(isi)
    if t > Tmax:
        ISIs_chunks.append(chunks.copy())
        chunks.clear()
        t = 0

ws = np.logspace(0, 2, 100)
spectrum_data = []
for w in ws:
    fws_real = []
    fws_img = []
    for isis in ISIs_chunks:
        fw = fourier_transformation_isis(w, isis)
        fws_real.append(fw.real)
        fws_img.append(fw.imag)
    spectrum_data.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))

mean = np.mean(ISIs_flat)
cv2 = np.var(ISIs_flat)/np.power(mean, 2)
print(np.sqrt(cv2), cv2)
Ps = []
ts = np.linspace(0, 5, 5000)
for t in ts:
    if t == 0:
        x = 0
    else:
        t3 = np.power(t, 3)
        x = np.sqrt(mean / (2 * np.pi * cv2 * t3)) * np.exp(-np.power(t - mean, 2) / (2 * mean * cv2 * t))
    Ps.append(x)

Ps_fourier = []
Fs = []
dt = ts[1] - ts[0]
for w in ws:
    P_fourier = 0
    for t, P in zip(ts ,Ps):
        P_fourier += P*np.exp(1j*w*t)*dt
    Fs.append(w)
    Ps_fourier.append(P_fourier)

spectrum_theory = []
for P_fourier in Ps_fourier:
    real = P_fourier.real
    imag = P_fourier.imag
    spectrum_f = 1/mean * (1 - (real**2 + imag**2))/((1-real)**2 + imag**2)
    spectrum_theory.append(spectrum_f)

rho1_mean = np.mean(correlations)
rho1_std = np.std(correlations)
rho1_std_error = rho1_std/np.sqrt(len(correlations))
print(np.mean(correlations), np.std(correlations))

ax.set_xlabel("$\omega$")
ax.set_ylabel("$S(\omega)$")
ax.set_xscale("log")
ax.set_yscale("log")

ax.plot(ws, spectrum_theory, c="k")
ax.plot(ws, spectrum_data, label=rf"$\rho_1 \pm \sigma_{{\rho_1}}/\sqrt{{N}}= {rho1_mean:.2f} \pm {rho1_std_error:.2f}$")
ax.legend()

plt.savefig(home + "/Desktop/Ca data/Spikes/HEK/Plots2/HEK2_bapta_power_spectrum.png")
plt.show()
