import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st
import functions as fc
import default_parameters as df

if __name__ == "__main__":
    cas = [0.50]
    cR = 0.2
    cT = 0.5
    for ci in cas:
        tau = 1.0
        p = 1.0
        K = 10
        M = 3
        N = 5

        """
        1. Load Data for the puff current from up to k_max=10 simulations and plot the variance
        with respect to the interval duration over which the current is averaged
        """
        home = os.path.expanduser("~")
        k_max = 1
        factor_max = 1000
        varss = []
        for _ in range(factor_max):
            vars = []
            for _ in range(k_max):
                vars.append([])
            varss.append(vars)

        dts = [None]*factor_max
        for k in range(k_max):
            data = df.load_traces_fixed_ci_markov(tau, p, ci)
            ts, cas, jpuffs, adap = np.transpose(data)

            print("mean sim: ", np.mean(jpuffs))
            dt0 = 0.01
            factors = np.arange(1, 1001)
            for idx, f in enumerate(factors):
                cg_list = fc.coarse_grain_list(jpuffs, f)
                dt = dt0 * f
                dts[idx] = dt
                varss[idx][k] = np.var(cg_list)*dt

        vars_mean_and_var = []
        for vars in varss:
            vars_mean_and_var.append([np.mean(vars), np.std(vars)])

        st.set_default_plot_style()
        fig = plt.figure(tight_layout=True, figsize=(4, 3))
        gs = gridspec.GridSpec(1, 1)
        ax = fig.add_subplot(gs[0, 0])
        st.remove_top_right_axis([ax])
        ax.plot(dts, [np.sqrt(x[0]) for x in vars_mean_and_var], label=r"$\sigma \sqrt{\Delta t}  =  \left\langle \left(\Delta \int_0^{\Delta t} dt\,  j_{\rm puff}(t)\right)^2 \right\rangle^{1/2}$ ")

        ax.set_ylabel("$\sigma \sqrt{\Delta t}$")
        ax.set_xlabel(r"$\Delta t$")
        ax.set_xscale("log")

        """
        2. Calculate the noise intensity from the theory
        """

        r_opn = df.r_opn(ci, s=1., N=N)
        r_ref = df.r_ref
        r_cls = df.r_cls

        D_single_theory = fc.noise_intensity_jp_single_theory(ci, N, M, s=1.0)

        ax.axhline(np.sqrt(2*K*D_single_theory), ls=":", c="C7", label="$\sqrt{{2 K D_x}}={:.2f}$".format(np.sqrt(2*K*D_single_theory)))
        ax.set_xlim([0.01, 10])
        legend = ax.legend(fancybox=False, framealpha=1.0, fontsize=9)
        ax.set_ylim([0., 1.])

        plt.show()
        plt.close()