import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from functions import *

if __name__ == "__main__":
    set_default_plot_style()
    cas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for ca in cas:
        ca_fix = ca
        ca_rest = 0.33
        N = 10
        m = 4
        n = 4

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
            folder = "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/"
            file =  "ca_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N{:d}_0.dat".format(ca_fix, N)
            data_spike = np.loadtxt(home + folder + file)

            ts, cas, jpuffs, adap = np.transpose(data_spike)
            print("mean sim: ", np.mean(jpuffs))
            dt0 = 0.1
            factors = np.arange(1, 101)
            for idx, f in enumerate(factors):
                cg_list = coarse_grain_list(jpuffs, f)
                dt = dt0 * f
                dts[idx] = dt
                varss[idx][k] = np.var(cg_list)*dt

        vars_mean_and_var = []
        for vars in varss:
            vars_mean_and_var.append([np.mean(vars), np.std(vars)])

        fig = plt.figure(tight_layout=True, figsize=(4, 3))
        gs = gridspec.GridSpec(1, 1)
        ax = fig.add_subplot(gs[0, 0])
        remove_top_right_axis([ax])
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

        r_opn = 0.13 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
        r_ref = 1.3 * np.power(ca_fix / ca_rest, 3) * ((1 + ca_rest**3)/(1 + ca_fix**3))
        r_cls = 50

        p0s = steady_states_theory_invert_A(r_ref, r_opn, r_cls, m, n)

        xs = [0, 0, 0, 0, 4, 3, 2, 1]
        idxs = [0, 1, 2, 3, 4, 5, 6, 7]
        mean = sum([x*p for x, p in zip(xs, p0s)])

        D_theory = 0
        for k in idxs:
            sum_over_i = 0
            f_from_k_to = f_from_k_invert_A(k, r_ref, r_opn, r_cls, m, n)
            for i in idxs:
                sum_over_i += xs[i]*f_from_k_to[i]
            D_theory += xs[k] * p0s[k]*sum_over_i

        ax.axhline(np.sqrt(2*N*D_theory), ls="--", c="C7", label="$\sqrt{{2D_N}}={:.2f}$".format(np.sqrt(2*N*D_theory)))
        ax.set_xlim([0.1, 10])
        legend = ax.legend(loc=3, fancybox=False, edgecolor="k", framealpha=1.0)
        legend.get_frame().set_linewidth(0.5)
        plt.savefig(home + "/Data/Calcium/Plots/4_markov_jpuff_D_var_over_deltat_static.pdf", transparent=True)
        plt.show()
        plt.close()