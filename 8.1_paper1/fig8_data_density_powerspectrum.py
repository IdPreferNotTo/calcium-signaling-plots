import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import os

import scipy.stats

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
    file_str = home + "/Data/calcium_spikes_experimental/Spikes/HEK/HEK2_bapta_ratio.dat"
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
    fig = plt.figure(tight_layout=True, figsize=(9, 2.5))
    gs = gridspec.GridSpec(1, 3)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    st.remove_top_right_axis([ax1, ax2, ax3])
    ax1.text(0.10, 0.95, "A", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.15, 0.95, "B", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.15, 0.95, "C", fontsize=13, transform=ax3.transAxes, va='top')

    ax1.set_xlim([0, 2])
    ax1.set_xlabel(r"$\tilde{T}$")
    ax1.set_ylabel(r"$p_{\rm ISI}(\tilde{T})$")

    mean = np.mean(ISIs_flat)
    std = np.std(ISIs_flat)
    var = np.power(std, 2)
    cv = std / mean
    cv2 = np.power(cv, 2)
    skew = scipy.stats.skew(ISIs_flat)
    print(skew)
    #print(n, mean, std, cv)
    tmax = max(ISIs_flat)
    ax1.hist(ISIs_flat, color=st.colors[4], bins=25, density=True)
    ts_inv_gau = np.linspace(0, 3, 1001)
    inv_gaus = get_inverse_gaussian(ts_inv_gau, mean, cv)
    ax1.plot(ts_inv_gau, inv_gaus, lw=1, c="k")


    fano_factors = []
    maxT = sum(ISIs_flat)
    spike_times = []
    spike_time = 0
    for isi in ISIs_flat:
        spike_time += isi
        spike_times.append(spike_time)

    means_count = []
    vars_count = []
    fanos_count = []
    dts = np.logspace(-1, 1, 1000)
    for dt in dts:
        mean_count, var_count = fc.fano_factor_interspike_intervals(ISIs_flat, dt)
        means_count.append(mean_count)
        vars_count.append(var_count)
        fanos_count.append(var_count/mean_count)

    ax2.plot(dts, fanos_count, c=st.colors[4], label=rf"$CV_T = {cv:.2f}$")
    ax2.axvline(mean, ls=":", c="C7")
    ax2.set_xscale("log")
    ax2.set_xlabel("$t$ / s")
    ax2.set_ylabel("$F(t)$")
    ax2.axhline(cv2, ls=":", c="C7")
    ax2.text(0.2, 0.08, "$CV_T^2$")
    ax2.set_ylim([0.0, 1.0])
    ax2.set_xlim([0.1, 10])
    ax2.legend(fancybox=False, loc=1, framealpha=1.)

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

    fs = np.logspace(-1, 1, 100)
    spectrum_data = []
    for f in fs:
        fws_real = []
        fws_img = []
        for isis in ISIs_chunks:
            fw = fc.fourier_transformation_isis(f, isis)
            fws_real.append(fw.real)
            fws_img.append(fw.imag)
        spectrum_data.append((1. / Tmax) * (np.var(fws_real) + np.var(fws_img)))

    spectrum_data = fc.power_spectrum_isis(fs, ISIs_flat, Tmax=25)
    Ps_fourier = []
    Fs = []
    dt = ts_inv_gau[1] - ts_inv_gau[0]
    for f in fs:
        P_fourier = 0
        for t, P in zip(ts_inv_gau, inv_gaus):
            P_fourier += P*np.exp(1j*2*np.pi*f*t)*dt
        Fs.append(fs)
        Ps_fourier.append(P_fourier)
    spectrum_theory = []
    for P_fourier in Ps_fourier:
        real = P_fourier.real
        imag = P_fourier.imag
        spectrum_f = 1/mean * (1 - (real**2 + imag**2))/((1-real)**2 + imag**2)
        spectrum_theory.append(spectrum_f)


    ax3.set_xlim([0.1, 10])
    ax3.set_xticks([0.1, 1, 10])
    ax3.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax3.axhline(1, ls=":", c="C7")
    ax3.axvline(1, ls=":", c="C7")
    ax3.axhline(cv**2, ls=":", c="C7")
    ax3.text(0.4, 1.3, "$r_0$")
    ax3.text(2, 1.3*cv**2, "$r_0 CV_T^2$")
    ax3.plot(fs, spectrum_theory, lw=1, c="k")
    ax3.plot(fs, spectrum_data, c=st.colors[4])
    ax3.set_xlabel("$f$ / s$^{-1}$")
    ax3.set_ylabel("$S(f)$")
    ax3.set_xscale("log")
    ax3.set_yscale("log")

    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig8.pdf", transparent=True)
    plt.show()
