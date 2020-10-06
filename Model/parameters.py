import json
import numpy as np
import os


def write_parameter_file(N, tau, ca_rest, ca_max, ip3, n_channel, ca_current, n_cluster, ipi_rest):
    home = os.path.expanduser('~')
    with open(home + '/Parameter/{}.json'.format(N), 'w') as file:
        json.dump({
            "num parameter": {
                "dt": 10 ** (-4),
                "spikes": 1_000,
            },
            "cell": {
                "timeconstant": tau,
                "calcium_basal": ca_rest,
                "calcium_spike": ca_max,
                "ip3": ip3,
            },
            "cluster": {
                "channel": {
                    "number": n_channel,
                    "current": ca_current,
                },
                "number": n_cluster,
                "interpuff interval": ipi_rest,
                "max_open_probability": 0.8,
                "actiavtion": {
                    "halfmax": 0.21,
                    "hill": 1.9,
                },
                "inhibition": {
                    "halfmax": 52 * np.power(1 + np.power(0.05 / ip3, 4), -1),
                    "hill": 3.9,
                },
            },
        },
            file, indent=0)


if __name__ == "__main__":
    N = 0
    for i in range(4):
        N += 1
        tau = 2 # Time constant of the Ca2+ removal into the ER and extracellular medium
        ca_rest = 0.05 # Basal Ca2+ concentration
        ca_max = 1 # Ca2+ concentration at which a spike a guaranteed to be triggered
        ip3 = 0.1

        n_channel = 5
        ca_current = 0.001
        n_cluster = 100
        ipi_rest = 1
        write_parameter_file(N, tau, ca_rest, ca_max, ip3, n_channel, ca_current, n_cluster, ipi_rest)

