import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import functions as fc

home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Spikes/HEK/HEK2_bapta_ratio.dat"
data = np.loadtxt(file_str)

ISIs = []
MEANs = []
STDs = []
n = len(data[0])
non_stat = [7, 10, 17, 19, 21, 22, 25, 26]
for j in range(1, n):
    ts = [x[0] for x in data]
    cas = [x[j] for x in data]

    spiking: bool = False
    t_tmp: float = 0
    spike_times = []
    if j in non_stat:
        continue
    elif j==2:
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
    ISIs.append(isis)
    n = len(isis)
    mean = np.mean(isis)
    MEANs.append(mean)
    std = np.std(isis)
    STDs.append(std)
    cv = std / mean
    cv2 = np.power(cv, 2)
    print(j, n, mean, std, cv)

fc.set_default_plot_style()
fig = plt.figure(tight_layout=True, figsize=(4, (3 / 2) * (6 / 2)))
gs = gridspec.GridSpec(3, 1)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])
fc.remove_top_right_axis([ax1, ax2, ax3])

ax1.set_xlabel(r"$\langle T_k \rangle _i$")
ax1.set_ylabel(r"$\sqrt{\langle \Delta T_k^2 \rangle _i}$")
ax1.scatter(MEANs, STDs, fc="w", ec="C0", s=20, label="HEK")
coef = np.polyfit(MEANs, STDs,1)
t_ref = -coef[1]/coef[0]
print(t_ref)

xs = np.linspace(0, 400, 100)
ys = coef[1] + coef[0]*xs
ax1.plot(xs, ys, ls=":", lw=1, c="k", label="Lin. fit")
ax1.set_ylim([0, 100])
ax1.set_xlim([0, 400])
legend = ax1.legend(loc=2, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 8})
legend.get_frame().set_linewidth(0.5)

ax2.set_xlim([0, 2])
ax2.set_xlabel(r"$T_{k, i} / \langle T_k \rangle _i$")
ax2.set_ylabel("$P(T_{k, i})$")
ISIs_no_ref = []
for isis in ISIs:
    isis_no_ref = []
    for isi in isis:
        isis_no_ref.append(isi - t_ref)
    ISIs_no_ref.append([t/np.mean(isis_no_ref) for t in isis_no_ref])

ISIs_flat = [item for sublist in ISIs_no_ref for item in sublist]

mean = np.mean(ISIs_flat)
std = np.std(ISIs_flat)
var = np.power(std, 2)
cv = std / mean
cv2 = np.power(cv, 2)
print(n, mean, std, cv)
tmax = max(ISIs_flat)
ax2.hist(ISIs_flat, bins=25, density=True, label="HEK")
ts = np.linspace(0.1, 2, 100)
Ps = []

for t in ts:
    if t == 0:
        x = 0
    else:
        t3 = np.power(t, 3)
        x = np.sqrt(mean / (2 * np.pi * cv2 * t3)) * np.exp(-np.power(t - mean, 2) / (2 * mean * cv2 * t))
    Ps.append(x)
ax2.plot(ts, Ps, lw=1, ls="--", c="k", label="Inv. Gaus")
legend = ax2.legend(loc=2, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 8})
legend.get_frame().set_linewidth(0.5)

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
        fw = fc.fourier_transformation_isis(w, isis)
        fws_real.append(fw.real)
        fws_img.append(fw.imag)
    spectrum_data.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))

Ps_fourier = []
Fs = []
dt = ts[1] - ts[0]
for w in ws:
    P_fourier = 0
    for t, P in zip(ts, Ps):
        P_fourier += P*np.exp(1j*w*t)*dt
    Fs.append(w)
    Ps_fourier.append(P_fourier)

spectrum_theory = []
for P_fourier in Ps_fourier:
    real = P_fourier.real
    imag = P_fourier.imag
    spectrum_f = 1/mean * (1 - (real**2 + imag**2))/((1-real)**2 + imag**2)
    spectrum_theory.append(spectrum_f)

ax3.plot(ws, spectrum_theory, lw=1, ls="--", c="k", label="Inv. Gaus")
ax3.plot(ws, spectrum_data, lw=1, label="HEK")
ax3.set_xlabel("$\omega$")
ax3.set_ylabel("$S(\omega)$")
ax3.set_xscale("log")
ax3.set_yscale("log")
legend = ax3.legend(loc=2, fancybox=False, edgecolor="k", framealpha=1.0, prop={'size': 8})
legend.get_frame().set_linewidth(0.5)

plt.savefig(home + f"/Paper/RamLin21_Ca/Plots/fig2_HEK_ISI_density_ps.pdf")
plt.show()
