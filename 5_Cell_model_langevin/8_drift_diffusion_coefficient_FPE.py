import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import os

import functions as fc
import styles as st


def mean_puff_singke_local(x, n, m, IP3, r_opn_single, r_ref):
    r_opn = n * r_opn_single * np.power(x / 0.2, 3) * ((1 + 0.2 ** 3) / (1 + x ** 3)) * np.power(IP3 / 1., 3) * ((1. + 1. ** 3) / (1. + IP3 ** 3))
    r_ref = r_ref
    r_cls = 50

    p0s = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, n, m)
    xs = fc.get_states(n, m)
    mean = sum([x * p for x, p in zip(xs, p0s)])
    return mean

if __name__ == "__main__":
    st.set_default_plot_style()
    fig = plt.figure(tight_layout=True, figsize=(6, 2.5))
    gs = gs.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    st.remove_top_right_axis([ax1, ax2])
    ax1.axhline(0, lw=1, color="C7")

    rcls = 50
    rref = 20
    ropn = 0.1
    N = 5
    K = 10
    M = 3
    cr = 0.200
    ct = 0.500
    cas = np.linspace(cr, ct, 100)
    taus = [5.62, 1.78]
    ps = [0.0126, 0.0355]
    for i, (tau, p) in enumerate(zip(taus, ps)):
        means = []
        jleaks = []
        jpuffs = []
        intensities = []
        for ca in cas:
            jleak =  -(ca-cr)/tau

            jpuff = K*p*mean_puff_singke_local(ca, N, M, 1, ropn, rref)
            mean = jleak + jpuff
            intensity = np.sqrt(2*K*p*p*fc.intensity_puff_single(ca, N, M, 1, ropn, rref))
            jleaks.append(jleak)
            jpuffs.append(jpuff)
            means.append(mean)
            intensities.append(intensity)

        ax1.plot(cas, jleaks, ls="--", c="k")
        ax1.plot(cas, jpuffs, ls=":", c=st.colors[i+1])
        ax1.plot(cas, means, c=st.colors[i+1])
        ax2.plot(cas, intensities, c=st.colors[i+1])
        ax1.axvline(0.5, c="k")
        ax1.axvline(0.5, c="k")
    home = os.path.expanduser("~")
    #plt.savefig(home + "/Data/Calcium/Plots/8_langevin_drift_diffusion_coef.pdf", transparent=True)
    plt.show()
