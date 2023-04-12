import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

import styles as st
import functions as fc
import default_parameters as df

def interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    return xc,yc


def exponential_cer(t, cer0, cer8, tau):
    return cer8 + (cer0 - cer8) * np.exp(-t / tau)


def turn_2d_isis_into_spike_times(ISIss):
    spike_times = np.zeros((6000, 25))
    for i, ISIs in enumerate(ISIss):
        t = 0
        for ii, I in enumerate(ISIs):
            t += I
            spike_times[i, ii] = t
    return spike_times


if __name__ == "__main__":
    # Set Plot style
    home = os.path.expanduser("~")
    st.set_default_plot_style()
    w = 3.25 * 1.25
    h = 2.25 * 1.25
    fig = plt.figure(tight_layout=True, figsize=(w, h))
    gs = gridspec.GridSpec(3, 1, height_ratios= [0.2, 0.4, 0.4], hspace=0)
    ax0 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[2])
    ax1 = fig.add_subplot(gs[1], sharex=ax2)
    st.remove_everything([ax0])
    st.remove_top_right_axis([ax1, ax2])

    ax1.set_ylabel("$r(t)$")
    ax1.set_xticklabels([])

    ax2.set_ylim([0.825, 1.05])
    ax2.set_ylabel(r"$\langle c_{\rm er}(t) \rangle$")
    ax2.set_xlabel(r"$t$")

    # Parameters
    tau = 5
    p = 0.015
    eps_er = 0.05
    tau_er = 200

    # Plot Stimulation
    ax0.set_ylim([-0.5, 1.5])
    ax0.plot([-500, 0, 0, 1000, 1000, 2000], [0, 0, 1, 1, 0, 0], c="k")
    ax0.text(-250, 1.07, "Stim.\ off", ha="center", va='bottom')
    ax0.text(500, 1.07, "Stim.\ on", ha="center", va='bottom')
    ax0.text(1500, 1.07, "Stim.\ off", ha="center", va='bottom')

    # Plot spike times by vertical lines
    data = df.load_traces_markov(tau, p, True, tau_er, eps_er)
    ts, cs, js, cers = np.transpose(data)
    ax2.plot(ts[:100_000], cers[:100_000], c="C7", lw=1, zorder=1)
    ts_spike = []
    for t, c, j in zip(ts, cs, js):
        if c == 0.5 and j == 0:
            ts_spike.append(t)
    for t_spike in ts_spike[:16]:
        ax0.plot([t_spike, t_spike], [0.2, 0.8], lw=1, c="k")

    # Plot r(t)
    data_isi_trans = df.load_spike_times_markov_transient(tau, p, tau_er, eps_er)
    spike_times_trans = turn_2d_isis_into_spike_times(data_isi_trans)
    spike_times_trans_flat = [item for sublist in spike_times_trans for item in sublist]
    ts = np.linspace(0, 1000, 100)
    n_over_t = np.zeros(100)
    dt = 1.
    for i, t in enumerate(ts):
        for ti in spike_times_trans_flat:
            if t-dt/2 < ti and ti < t+dt/2:
                n_over_t[i] += 1
    r_over_t = n_over_t / (dt * 6_000)
    ts = np.append(ts, [1_000, 2_000])
    r_over_t = np.append(r_over_t, [0, 0])
    ax1.plot([-500, 0], [0, 0], lw=1, zorder=2, c="C7")
    ax1.plot(ts, r_over_t, lw=1, zorder=2, c="C7")
    ax1.axvline(0, ls=":", c="C7")
    ax1.axvline(1000, ls=":", c="C7")

    # Plot c_er(t)
    file = f"transient_adaptation_markov_ip1.00_taua{tau_er:.2e}_ampa{eps_er:.2e}_tau{tau:.2e}_j{p:.2e}_K10_5.dat"
    data = np.loadtxt("/home/lukas/Data/calcium_spikes_markov/Data_adap_transient/" + file,
                      usecols=np.arange(0, 1000))
    mean_cers_on = np.mean(data, axis=0)
    ts = np.linspace(0, 1000, 1000)
    popt_on, pcov = curve_fit(fc.exponential_cer, ts, mean_cers_on, p0=(0.9, 100))
    mean_cers_on_fit = popt_on[0] + (1 - popt_on[0]) * np.exp(-ts / popt_on[1])

    mean_cer_off1 = mean_cers_on_fit[-1]
    mean_cer_off = mean_cer_off1
    mean_cers_off = []
    dt = ts[1] - ts[0]
    for t in ts:
        mean_cers_off.append(mean_cer_off)
        dcer = -(mean_cer_off - 1.) / tau_er
        mean_cer_off += dcer*dt
    popt_off, pcov = curve_fit(exponential_cer, ts, mean_cers_off, p0=(0.9, 1.0, 100))
    mean_cers_off_fit = popt_off[1] + (popt_off[0] - popt_off[1]) * np.exp(-ts / popt_off[2])

    ax2.annotate(r'$\tau_{\rm eff}$', xy=(0, 0.95),
                 xytext=(-250, 0.9), va="top", ha="right",
                 arrowprops=dict( arrowstyle="->" ))
    ax2.annotate(r'$\tau_{\rm er}$', xy=(1250, 0.95),
                 xytext=(1500, 0.9), va="top", ha="left",
                 arrowprops=dict( arrowstyle="->" ))

    ax2.plot(ts, mean_cers_on_fit, c="k", ls="--", label="Fit " + rf"$\tau_{{\rm eff}} = {popt_on[1]:.0f}$", zorder=3)
    ax2.axvline(0, ls=":", c="C7")
    ax2.axvline(1000, ls=":", c="C7")

    ts = np.insert(ts, 0, -500)
    mean_cers_on = np.insert(mean_cers_on, 0, 1)
    mean_cers_off = np.insert(mean_cers_off, 0, mean_cer_off1)
    ax2.plot(ts, mean_cers_on, color=st.colors[1])
    ax2.plot(ts + 1000, mean_cers_off, color=st.colors[1])

    data_rate_cer = np.loadtxt(home + f"/Data/calcium_spikes_theory/cer_r0_langevin_tau{tau:.2e}_j{p:.2e}.dat")
    cer_fix, r_fix = np.transpose(data_rate_cer)
    rates = []
    for cer in mean_cers_on[1:]:
        r = fc.linear_interpolate(cer, cer_fix, r_fix)
        rates.append(r)

    ts = np.insert(ts, 1, 0) # very dirty ...
    ts = np.append(ts, [1000, 2000])
    rates = np.insert(rates, 0, [0, 0])
    rates = np.append(rates, [0, 0])
    ax1.plot(ts, rates, color=st.colors[1], zorder = 3)
    plt.savefig(home + f"/Dropbox/LUKAS_BENJAMIN/RamLin22_2_BiophysJ/figures/fig5.pdf")
    plt.show()

