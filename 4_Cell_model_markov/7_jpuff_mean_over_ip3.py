import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import functions as fc

if __name__ == "__main__":

    cas_left = np.linspace(0.23, 0.33, 10)
    cas = np.linspace(0.33, 1, 100)
    cas_right = np.linspace(1, 1.1, 10)
    jpuff_over_ca_left = []
    jpuff_over_ca = []
    jpuff_over_ca_right = []
    caR = 0.33
    j = 0.1
    N = 10
    n = 4
    m = 3
    rcls = 50
    ropn_max = 0.13 * (1 + np.power(caR, 3))/np.power(caR, 3)
    rref_max = 1.3 * (1 + np.power(caR, 3))/np.power(caR, 3)
    alpha = 3
    beta = 3
    def ropn(x, y=1):
        return n*ropn_max*(np.power(x, alpha)/(1. + np.power(x, alpha)))*(np.power(y, beta)/(1. + np.power(y, beta)))

    def rref(x, y=1):
        return n*rref_max*(np.power(x, alpha)/(1. + np.power(x, alpha)))*(np.power(y, beta)/(1. + np.power(y, beta)))

    def t_opn():
        return (n+1)/(2*rcls)

    def t_cls(x, y=1):
        return 1/ropn(x, y) + m/rref(x, y)

    for ca in cas_left:
        jpuff = j*N*((n+2)/3)*(t_opn() / (t_opn() + t_cls(ca)))
        jpuff_over_ca_left.append(jpuff)
    for ca in cas:
        jpuff = j*N*((n+2)/3)*(t_opn() / (t_opn() + t_cls(ca)))
        jpuff_over_ca.append(jpuff)
    for ca in cas_right:
        jpuff = j*N*((n+2)/3)*(t_opn() / (t_opn() + t_cls(ca)))
        jpuff_over_ca_right.append(jpuff)

    fc.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(4, 6 / 2))
    grids = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(grids[:])
    axins = ax.inset_axes((0.6, 0.2, 0.4, 0.4))
    fc.remove_top_right_axis([ax])

    axins.plot(cas_left, jpuff_over_ca_left, ls=":", c="C0")
    axins.plot(cas, jpuff_over_ca, c="C0")
    axins.plot(cas_right, jpuff_over_ca_right, ls=":", c="C0")
    axins.axvline(0.33, ls=":", c="C7")
    axins.axvline(1.00, ls=":", c="C7")
    axins.set_xlim([0.23, 1.1])
    axins.set_xticks([0.33, 1])
    axins.set_xticklabels(["$c_R$", "$c_T$"])
    axins.set_ylim([0, 0.3])
    axins.set_xlabel("c")
    axins.set_ylabel(r"$j_{\rm puff}(c)$")

    dca = cas[1] - cas[0]
    jpuff_mean = 1/(1 - 0.33) * sum(jpuff_over_ca)*dca

    axins.plot((0.33, 1), (jpuff_mean, jpuff_mean), ls="--", c="C0")

    jpuff_mean_over_ip3 = []
    jpuff_mean_over_ip3_at_caR = []
    ip3s = np.linspace(0, 2, 100)
    for ip3 in ip3s:
        jpuff_over_ca = []
        for ca in cas:
            jpuff = j * N * ((n + 2) / 3) * (t_opn() / (t_opn() + t_cls(ca, ip3)))
            if ca == 0.33:
                jpuff_mean_over_ip3_at_caR.append(jpuff)
            jpuff_over_ca.append(jpuff)
        jpuff_mean = 1 / (1 - 0.33) * sum(jpuff_over_ca) * dca
        jpuff_mean_over_ip3.append(jpuff_mean)

    ax.plot(ip3s, jpuff_mean_over_ip3_at_caR, c="C0", ls=":")
    ax.plot(ip3s, jpuff_mean_over_ip3)
    ax.set_xlabel(r"IP$_3/K_{\rm IP_3}$")
    ax.set_ylabel(r"$\overline{j_{\rm puff}}$ / [a.u.]")
    plt.show()