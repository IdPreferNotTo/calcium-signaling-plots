import json
import os


def print_parameter_files(N, M, tau_x, mu_x, D_x, delta_ca, tau_ca, alpha_ca, delta_inh, tau_inh, beta_inh, stepsize= 10**(-5), max_spikes=5_000):
    i = 1
    home = os.path.expanduser('~')
    with open(home + '/CLionProjects/calcium-phd/data/{}.json'.format(i), 'w') as file:
        json.dump({
            "numerical_parameter": {
                "stepsize": stepsize,
                "max_spikes": max_spikes,
            },
            "puff": {
                "sites":N,
                "channels":M,
                "tau":{
                    "mean": tau_x,
                    "stddev":0,
                },
                "mu":{
                    "mean": mu_x,
                    "stddev":0,
                },
                "D":{
                    "mean": D_x,
                    "stddev":0,
                },

            },
            "calcium": {
                "blip increment": delta_ca,
                "tau": tau_ca,
                "alpha": alpha_ca,
            },
            "inhibitor": {
                "spike increment": delta_inh,
                "tau": tau_inh,
                "beta": beta_inh,
            },
        },
            file, indent=0)


if __name__ == "__main__":
    N = 15
    M = 6
    tau_x = 1
    mu_x = 1/6
    D_x = 0.5
    delta_ca = 0.01
    tau_ca = 4
    alpha_ca = 1
    delta_inh = 1.
    tau_inh = 10
    beta_inh = 0
    print_parameter_files(N, M, tau_x, mu_x, D_x, delta_ca, tau_ca, alpha_ca, delta_inh, tau_inh, beta_inh)