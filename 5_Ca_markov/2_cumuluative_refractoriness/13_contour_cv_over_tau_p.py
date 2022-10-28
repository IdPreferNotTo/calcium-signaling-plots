import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    pmax = 0.1
    pmin = 0.001
    tmax = 100
    tmin = 1
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(8, 3))
    gs = gridspec.GridSpec(1, 2)
    ax_cv = fig.add_subplot(gs[0])
    ax_dcv = fig.add_subplot(gs[1])
    ax_cv.text(0.1, 0.1, r"A", fontsize=13, transform=ax_cv.transAxes)
    ax_dcv.text(0.1, 0.1, r"B", fontsize=13, transform=ax_dcv.transAxes)

    ax_cv.text(0.65, 0.65, "mean-driven", fontsize=13, rotation=-41, transform=ax_dcv.transAxes, ha="center", va="center")
    ax_cv.text(0.25, 0.25, "excitable", fontsize=13, rotation=-41, transform=ax_dcv.transAxes, ha="center", va="center")
    ax_cv.arrow(0.33, 0.33, 0.1, 0.1 * np.tan((50 / 360) * 2 * np.pi), transform=ax_dcv.transAxes, fc="k", length_includes_head=True, head_width=0.05 / 2, head_length=0.1 / 2, )
    colors = st.Colors()

    size = 41
    K = 10
    N = 5
    M = 3
    home = os.path.expanduser("~")
    data_renewal = np.loadtxt(home + f"/Data/calcium_spikes_theory/markov_isi_mean_CV_K{K:d}_N{N:d}_no_adap.dat")
    taus, js, Ts_renew, cvs_renew, ns = np.transpose(data_renewal)

    taus = taus.reshape(size, size)
    js = js.reshape(size, size)
    Ts_renew = Ts_renew.reshape(size, size)
    cvs_renew = cvs_renew.reshape(size, size)
    cmap_cividis = plt.get_cmap("YlGnBu", 10)
    contour_cv = ax_cv.pcolormesh(taus, js, cvs_renew, linewidth=0, rasterized=True, shading='gouraud', vmin=0., vmax=1.0, cmap=cmap_cividis)
    divider = make_axes_locatable(ax_cv)
    cax_cv = divider.append_axes('right', size='5%', pad=0.05)
    cbar_cv = fig.colorbar(contour_cv, cax=cax_cv, orientation='vertical')
    cbar_cv.set_label(r"$\langle T \rangle$ / s", loc="center")
    ax_cv.set_ylabel(r"$p$")
    ax_cv.set_xlabel(r"$\tau$")
    ax_cv.set_xscale("log")
    ax_cv.set_yscale("log")
    ax_cv.set_xlim([tmin, tmax])
    ax_cv.set_ylim([pmin, pmax])

    tau_er = 200
    eps_er = 0.1
    data_cumulative_ref = np.loadtxt(home + f"/Data/calcium_spikes_theory/markov_ca_stationary_statistics_tau_er{tau_er:.2e}_eps_er{eps_er:.2e}_K{K:d}_N{N:d}_adap.dat")
    taus, js, Ts_refrac, cvs_refrac, ns = np.transpose(data_renewal)
    taus = taus.reshape(size, size)
    js = js.reshape(size, size)
    Ts_refrac = Ts_refrac.reshape(size, size)
    cvs_refrac = cvs_refrac.reshape(size, size)
    dcvs = (cvs_refrac - cvs_renew)/(cvs_renew)
    ax_dcv.set_ylabel(r"$p$")
    ax_dcv.set_xlabel(r"$\tau$")
    ax_dcv.set_xscale("log")
    ax_dcv.set_yscale("log")
    ax_dcv.set_xlim([tmin, tmax])
    ax_dcv.set_ylim([pmin, pmax])

    cmap_coolwarm = plt.get_cmap("coolwarm", 10)
    cs_dcv_markov = ax_dcv.pcolormesh(taus, js, dcvs, linewidth=0, rasterized=True, shading='gouraud', vmin=-1.0,
                                      vmax=1.0, cmap=cmap_coolwarm)
    # ax_cv.contour(taus_adap, epss_adap, CVs, levels=[0])
    divider = make_axes_locatable(ax_dcv)
    cax_dcv_markov = divider.append_axes('right', size='5%', pad=0.05)
    cbar_dcv_markov = fig.colorbar(cs_dcv_markov, cax=cax_dcv_markov, orientation='vertical')
    cbar_dcv_markov.set_label(r"$(CV_T - CV_T^*)/CV_T^*$", loc="center")
    cbar_dcv_markov.set_ticks([-1., 0, 1.0])

    r_cls = 50
    r_ref = 20
    r_opn_single = 0.1
    ropnmax = 5 * r_opn_single * ((1. + np.power(0.20, 3)) / np.power(0.20, 3))
    r_opn_ct = ropnmax * (np.power(0.5, 3) / (1. + np.power(0.5, 3)))
    mean_puff = (6) * (7) / (6 * r_cls)
    tau_tot = 1 / r_opn_ct + 2 / r_ref + 6 / (2 * r_cls)
    mean_puff_ct = mean_puff / tau_tot
    ps_bif = []
    ps_bif2 = []
    ps_bif3 = []
    for tau in taus:
        c = (0.5 - 0.2) / (10 * mean_puff_ct * tau)
        ps_bif.append(1.0*c)
        ps_bif2.append((3/4)*c)
        ps_bif3.append((4/3)*c)

    ax_cv.plot(taus, ps_bif, ls="--", c=colors.palette[5])
    ax_dcv.plot(taus, ps_bif, ls="--", c=colors.palette[5])
    plt.show()


