import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import styles as st
import functions as fc

if __name__ == "__main__":
    cas = [0.33] #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for ca in cas:
        ca_fix = 0.99
        ca_rest = 0.33
        tau = 1.0
        j = 1.0
        K = 10
        M = 3
        N = 5

        """
        1. Load Data for the puff current from up to k_max=10 simulations and plot the variance
        with respect to the interval duration over which the current is averaged
        """
        home = os.path.expanduser("~")
        k_max = 1
        factor_max = 100
        varss = []
        for _ in range(factor_max):
            vars = []
            for _ in range(k_max):
                vars.append([])
            varss.append(vars)

        dts = [None]*factor_max
        for k in range(k_max):
            folder = "/Data/calcium_spikes_markov/ca_fix/"
            file =  f"ca_markov_cafix{ca_fix:.2f}_ip1.00_tau{tau:.2e}_j{j:.2e}_K{K:d}_{N:d}.dat"
            data_spike = np.loadtxt(home + folder + file)

            ts, cas, jpuffs, adap = np.transpose(data_spike)
            print("mean sim: ", np.mean(jpuffs))
            dt0 = 0.1
            factors = np.arange(1, 101)
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
        var_upper = [x[0] + 3*x[1] for x in vars_mean_and_var]
        var_lower = [x[0] - 3*x[1] for x in vars_mean_and_var]
        ax.fill_between(dts, var_lower, var_upper, alpha=0.2)
        #ax.set_ylim([0., 0.5])
        ax.set_ylabel("$\sigma \sqrt{\Delta t}$")
        ax.set_xlabel(r"$\Delta t$")
        ax.set_xscale("log")

        """
        2. Calculate the noise intensity from the theory
        """

        r_opn = N*0.2 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
        r_ref = 5 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
        r_cls = 50

        p0s = fc.steady_states_theory_invert_M(r_ref, r_opn, r_cls, N, M)
        xs = fc.get_states(N, M)
        mean = sum([x*p for x, p in zip(xs, p0s)])

        D_theory = 0
        for k in range(N+M):
            sum_over_i = 0
            f_from_k_to = fc.f_from_k_invert_M(k, r_ref, r_opn, r_cls, N, M)
            for i in range(N+M):
                sum_over_i += xs[i]*f_from_k_to[i]
            D_theory += xs[k] * p0s[k]*sum_over_i

        ax.axhline(np.sqrt(2*N*D_theory), ls="--", c="C7", label="$\sqrt{{2D_N}}={:.2f}$".format(np.sqrt(2*N*D_theory)))
        ax.set_xlim([0.1, 10])
        legend = ax.legend(loc=3, fancybox=False, edgecolor="k", framealpha=1.0)
        legend.get_frame().set_linewidth(0.5)
        #plt.savefig(home + "/Data/Calcium/Plots/puff_current_mean_var_static.pdf", transparent=True)
        plt.show()
        plt.close()