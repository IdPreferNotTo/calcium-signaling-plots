import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import os

import styles as st
import functions as fc


def get_inverse_gaussian(ts, mean, cv):
    inv_gaus = []
    for t in ts:
        if t == 0:
            p = 0
        else:
            p = np.sqrt(mean / (2 * np.pi * np.power(cv,2) * (t ** 3))) * np.exp(
                -(t - mean) ** 2 / (2 * mean * np.power(cv,2) * t))
        inv_gaus.append(p)
    return inv_gaus


if __name__ == "__main__":
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
        #print(j, n, mean, std, cv)
    ISIs_flat = [item/np.mean(sublist) for sublist in ISIs for item in sublist]

    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 2))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    st.remove_top_right_axis([ax1, ax2])

    ax1.set_xlim([0, 2])
    ax1.set_xlabel(r"$T_{i, k} / \langle T_i \rangle _k$")
    ax1.set_ylabel("$P(T_{i, k})$")

    mean = np.mean(ISIs_flat)
    std = np.std(ISIs_flat)
    var = np.power(std, 2)
    cv = std / mean
    cv2 = np.power(cv, 2)
    #print(n, mean, std, cv)
    tmax = max(ISIs_flat)
    ax1.hist(ISIs_flat, color=st.colors[4], bins=25, density=True, label="HEK")

    ts_inv_gau = np.linspace(0, 2, 1001)
    inv_gaus = get_inverse_gaussian(ts_inv_gau, mean, cv)
    ax1.plot(ts_inv_gau, inv_gaus, lw=1, c="k")


    ISIs_chunks = []
    chunks = []
    t = 0
    Tmax = 20
    for isi in ISIs_flat:
        t += isi
        chunks.append(isi)
        if t > Tmax:
            ISIs_chunks.append(chunks.copy())
            chunks.clear()
            t = 0

    covars = []
    len(ISIs_flat)
    var_ISI = np.var(ISIs_flat)
    for i in range(0, len(ISIs_flat), len(ISIs_flat)):
        I = ISIs_flat[i:i+len(ISIs_flat)]
        covar = fc.k_corr(I, I, 1)/var_ISI
        covars.append(covar)
    print(np.mean(covars), np.std(covars))

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

    spectrum_data = fc.power_spectrum_isis(ws, ISIs_flat, Tmax=25)
    Ps_fourier = []
    Fs = []
    dt = ts_inv_gau[1] - ts_inv_gau[0]
    for w in ws:
        P_fourier = 0
        for t, P in zip(ts_inv_gau, inv_gaus):
            P_fourier += P*np.exp(1j*w*t)*dt
        Fs.append(w)
        Ps_fourier.append(P_fourier)
    spectrum_theory = []
    for P_fourier in Ps_fourier:
        real = P_fourier.real
        imag = P_fourier.imag
        spectrum_f = 1/mean * (1 - (real**2 + imag**2))/((1-real)**2 + imag**2)
        spectrum_theory.append(spectrum_f)


    ax2.set_xlim([1, 100])
    ax2.set_xticks([1, 10, 100])
    ax2.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax2.axvline(2*np.pi, ls=":", c="C7")
    ax2.axhline(1, ls="--", c="C7")
    ax2.axhline(cv**2, ls="--", c="C7")
    ax2.text(2, 1.1, "$r_0$")
    ax2.text(20, 1.1*cv**2, "$r_0 C_V^2$")
    ax2.plot(ws, spectrum_theory, lw=1, c="k")
    ax2.plot(ws, spectrum_data, lw=1, c=st.colors[4])
    ax2.set_xlabel("$\omega$")
    ax2.set_ylabel("$S(\omega)$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig8.png", transparent=True)
    plt.show()
