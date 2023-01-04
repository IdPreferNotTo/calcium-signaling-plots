import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.ticker as ticker
import os
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    ieffs = []
    T0s = []
    T8s = []
    ISIs = []
    N = 0
    n = 0
    with open(home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/spike_times.dat", "r") as csvfile:
        sreader = csv.reader(csvfile, delimiter=' ')
        for line in sreader:
            ieffs.append(float(line[0]))
            T0s.append(float(line[1]))
            T8s.append(float(line[2]))
            isis = []
            for isi in line[3:]:
                isis.append(float(isi))
            imin = int(1.5*float(line[0]))
            n += 1
            N += len(isis)
            isi_no_transient = isis[imin+1:]
            ISIs.append(isi_no_transient)

    print(N, n)
    print(np.mean([int(1.5*x) for x in ieffs]))

    ISIs_flat1 = [item for sublist in ISIs for item in sublist]
    ISIs_flat = [item/np.mean(sublist) for sublist in ISIs for item in sublist]

    print(np.mean(ISIs_flat1), np.std(ISIs_flat1)/np.mean(ISIs_flat1))
    print(np.mean(ISIs_flat), np.std(ISIs_flat) / np.mean(ISIs_flat))
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(8, 5))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0 ,0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    st.remove_top_right_axis([ax2, ax3, ax4])
    ax1.text(0.10, 0.95, "A", fontsize=13, transform=ax1.transAxes, va='top')
    ax2.text(0.10, 0.95, "B", fontsize=13, transform=ax2.transAxes, va='top')
    ax3.text(0.15, 0.95, "C", fontsize=13, transform=ax3.transAxes, va='top')
    ax4.text(0.15, 0.95, "D", fontsize=13, transform=ax4.transAxes, va='top')

    size = 41
    N = 10
    n = 5
    m = 3
    home = os.path.expanduser("~")
    data_markov = np.loadtxt(home + f"/Data/calcium_spikes_theory/markov_ca_mean_CV_K{N:d}_N{n:d}_no_adap.dat")
    taus_m, js_m, Ts_m, cvs_m, num_m = np.transpose(data_markov)
    taus = np.logspace(0, 2, size)
    js = np.logspace(-3, -1, size)
    ISI_markov = np.empty([size, size])
    CV_markov = np.empty([size, size])
    for k, T in enumerate(Ts_m):
        if T >= 1_000:
            ISI_markov[k // size, k % size] = np.nan
        else:
            ISI_markov[k // size, k % size] = T
    for k, cv in enumerate(cvs_m):
        if cv == 1:
            CV_markov[k // size, k % size] = np.nan
        else:
            CV_markov[k // size, k % size] = cv

    colors = st.Colors()
    pmax = 0.03
    pmin = 0.003
    tmax = 50
    tmin = 5
    ax1.set_ylabel(r"$p$")
    ax1.set_xlabel(r"$\tau$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax1.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax1.set_xticklabels(["$5 \cdot 10^0$", "$5 \cdot 10^1$"])
    ax1.set_yticklabels(["$3 \cdot 10^{-3}$", "$3 \cdot 10^{-2}$"])
    ax1.set_xticks([1, 10])
    ax1.set_yticks([0.1, 1])
    ax1.set_xlim([1, 10])
    ax1.set_ylim([0.1, 1])
    cmap_cividis = plt.get_cmap("YlGnBu", 10)
    cs_cv_markov = ax1.pcolormesh([t/5 for t in taus], [j/0.03 for j in js], CV_markov, linewidth=0, rasterized=True, shading='gouraud', vmin=0., vmax=0.3, cmap=cmap_cividis)
    cs1 = ax1.contour([t/5 for t in taus], [j/0.03 for j in js], CV_markov, linewidths=1, levels=[0.15], colors="k")
    cs2 = ax1.contour([t/5 for t in taus], [j/0.03 for j in js], ISI_markov, linewidths=1, levels= [157], colors=st.colors[5])

    divider = make_axes_locatable(ax1)
    cax_cv_markov = divider.append_axes('right', size='5%', pad=0.05)
    cbar_cv_markov = fig.colorbar(cs_cv_markov, cax=cax_cv_markov, orientation='vertical')
    cbar_cv_markov.set_label(r"$CV_T$", loc="center")
    cbar_cv_markov.set_ticks([0, 0.15, 0.3])

    ax2.set_xlim([0, 2])
    ax2.set_ylim([0, 4])
    ax2.set_xlabel(r"$T$")
    ax2.set_ylabel(r"$p_{\rm ISI}(T)$")

    mean = np.mean(ISIs_flat)
    std = np.std(ISIs_flat)
    var = np.power(std, 2)
    cv = std / mean
    cv2 = np.power(cv, 2)

    tmax = max(ISIs_flat)
    ax2.hist(ISIs_flat, color=st.colors[5], bins=20, alpha=1.0, density=True, zorder=2, label="HEK cells")
    ts_inv_gau = np.linspace(0, 3, 1001)
    inv_gaus = get_inverse_gaussian(ts_inv_gau, mean, cv)
    #ax1.plot(ts_inv_gau, inv_gaus, lw=1, c="k")

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

    ax3.plot(dts, fanos_count, c=st.colors[5], label="Data")
    ax3.axvline(mean, ls=":", c="C7")
    ax3.set_xscale("log")
    ax3.set_xlabel("$t$")
    ax3.set_ylabel("$F(t)$")
    ax3.axhline(cv2, ls=":", c="C7")
    ax3.text(0.2, 0.08, "$CV_T^2$")
    ax3.set_ylim([0.0, 1.0])
    ax3.set_xlim([0.1, 10])

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

    spectrum_data = fc.power_spectrum_isis(fs, ISIs_flat, Tmax=10)
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


    ax4.set_xlim([0.1, 10])
    ax4.set_xticks([0.1, 1, 10])
    ax4.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax4.axhline(1, ls=":", c="C7")
    ax4.axvline(1, ls=":", c="C7")
    ax4.axhline(cv ** 2, ls=":", c="C7")
    ax4.text(0.4, 1.3, "$r_0$")
    ax4.text(2, 1.3 * cv ** 2, "$r_0 CV_T^2$")
    #ax3.plot(fs, spectrum_theory, lw=1, c="k")
    ax4.plot(fs, spectrum_data, c=st.colors[5])
    ax4.set_xlabel("$f$")
    ax4.set_ylabel("$S(f)$")
    ax4.set_xscale("log")
    ax4.set_yscale("log")

    tau = 14.446540880503141
    p = 0.004765556953179596
    ax1.scatter(tau/5, p/0.03, zorder=3, s=50, ec="k", fc=st.colors[5])
    folder = home + "/Data/calcium_spikes_markov/Data_no_adap_zoom/"
    file_spikes = f"spike_times_markov_ip1.00_tau{tau:.2e}_j{p:.2e}_K10_0.dat"
    isis_model = np.loadtxt(folder + file_spikes)
    mean_isi = np.mean(isis_model)
    cv_isi = np.std(isis_model)/mean_isi
    print(mean_isi, cv_isi)
    isis_model = [isi/mean_isi for isi in isis_model]
    fs = np.logspace(-1, 1, 100)
    spectrum_data = fc.power_spectrum_isis(fs, isis_model, Tmax=10)

    dts = np.logspace(-1, 1, 1000)
    fano = []
    for dt in dts:
        mean_count, var_count = fc.fano_factor_interspike_intervals(isis_model, dt)
        fano.append(var_count / mean_count)

    ax2.hist(isis_model, bins=20, alpha = 1.0, color="k", histtype=u'step', zorder=3, density=True, label="Model")
    ax3.plot(dts, fano, lw = 1, color="k", zorder=2, label="Model")
    ax4.plot(fs, spectrum_data, lw = 1, color="k")
    ax2.legend(fancybox=False, loc=1, fontsize=9, framealpha=1.)

    #plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_1_BiophysJ/figures/fig8.pdf", transparent=True)
    plt.show()
