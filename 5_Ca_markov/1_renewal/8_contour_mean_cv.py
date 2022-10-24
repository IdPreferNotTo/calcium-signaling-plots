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
    ax_mean_markov = fig.add_subplot(gs[0])
    ax_cv_markov = fig.add_subplot(gs[1])
    ax_mean_markov.text(0.1, 0.1, r"A", fontsize=13, transform=ax_mean_markov.transAxes)
    ax_cv_markov.text(0.1, 0.1, r"B", fontsize=13, transform=ax_cv_markov.transAxes)

    ax_cv_markov.text(0.65, 0.65, "mean-driven", fontsize=13, rotation=-41, transform=ax_cv_markov.transAxes, ha="center", va="center")
    ax_cv_markov.text(0.25, 0.25, "excitable", fontsize=13, rotation=-41, transform=ax_cv_markov.transAxes, ha="center", va="center")
    ax_cv_markov.arrow(0.33, 0.33, 0.1, 0.1*np.tan((50/360) * 2 * np.pi), transform=ax_cv_markov.transAxes, fc="k", length_includes_head=True,  head_width=0.05/2, head_length=0.1/2,)
    colors = st.Colors()

    size = 41
    K = 10
    N = 5
    M = 3
    home = os.path.expanduser("~")
    data_markov = np.loadtxt(home + f"/Data/calcium_spikes_theory/markov_ca_mean_CV_K{K:d}_N{N:d}_no_adap.dat")
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
    # st.remove_top_right_axis([ax11, ax12, ax13, ax14, ax21, ax22, ax23, ax24])

    cmap_cividis = plt.get_cmap("YlGnBu", 10)

    cs_mean_markov = ax_mean_markov.pcolormesh(taus, js, ISI_markov, linewidth=0, rasterized=True, shading='gouraud', cmap=cmap_cividis,
                                               norm=mcolors.SymLogNorm(linthresh=0.01, base=10, vmin=1., vmax=1000))
    #ax_cv_markov.contour(taus, js, ISI_markov, linewidths=1, levels= [157], colors=colors.palette[5])

    divider = make_axes_locatable(ax_mean_markov)
    cax_mean_markov = divider.append_axes('right', size='5%', pad=0.05)
    cbar_mean_markov = fig.colorbar(cs_mean_markov, cax=cax_mean_markov, orientation='vertical')
    cbar_mean_markov.set_label(r"$\langle T \rangle$ / s", loc="center")
    ax_mean_markov.set_ylabel(r"$p$")
    ax_mean_markov.set_xlabel(r"$\tau$")
    ax_mean_markov.set_xscale("log")
    ax_mean_markov.set_yscale("log")
    ax_mean_markov.set_xlim([tmin, tmax])
    ax_mean_markov.set_ylim([pmin, pmax])

    ax_cv_markov.set_ylabel(r"$p$")
    ax_cv_markov.set_xlabel(r"$\tau$")
    ax_cv_markov.set_xscale("log")
    ax_cv_markov.set_yscale("log")
    ax_cv_markov.set_xlim([tmin, tmax])
    ax_cv_markov.set_ylim([pmin, pmax])

    cs_cv_markov = ax_cv_markov.pcolormesh(taus, js, CV_markov, linewidth=0, rasterized=True, shading='gouraud', vmin=0., vmax=1.0, cmap=cmap_cividis)
    divider = make_axes_locatable(ax_cv_markov)
    cax_cv_markov = divider.append_axes('right', size='5%', pad=0.05)
    cbar_cv_markov = fig.colorbar(cs_cv_markov, cax=cax_cv_markov, orientation='vertical')
    cbar_cv_markov.set_label(r"$CV_T$", loc="center")
    cbar_cv_markov.set_ticks([0, 0.5, 1.0])
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

    ax_mean_markov.plot(taus, ps_bif, ls="--", c=colors.palette[5])
    ax_cv_markov.plot(taus, ps_bif, ls="--", c=colors.palette[5])
    plt.show()


