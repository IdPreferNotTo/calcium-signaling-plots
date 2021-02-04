import json
import numpy as np
import os


def write_parameter_file(N, i, bool_ca_fix, ca_fix, tau, current, bool_adap, tau_adap, amp_adap, ip3, number_cluster, mean_num_channel, num_closed_states, ratio, mean_ibi, r_close):
    home = os.path.expanduser('~')
    with open(home + '/CLionProjects/calcium-spikes-from-puff-phd/parameter/{}.json'.format(N), 'w') as file:
        json.dump({
            "num parameter": {
                "run": i,
                "output": "local",
                "dt": 10 ** (-4),
                "t_out": 10 ** (-1),
                "max puffs": 100_000,
                "max time": 10_000
            },
            "cell": {
                "timeconstant": tau,
                "blip current": current,
                "calcium rest": 0.33,
                "calcium threshold": 1,
                "ip3": ip3,
                "num cluster": number_cluster
            },
            "calcium fix":{
                "on":bool_ca_fix,
                "value":ca_fix,
            },
            "adaptation":{
                "on": bool_adap,
                "timeconstant": tau_adap,
                "amplitude": amp_adap,
            },
            "cluster": {
                "mean num channel": mean_num_channel,
                "number closed states": num_closed_states,
                "rate open single": (1 + ratio * (num_closed_states-1)) / mean_ibi,
                "rate ref single": (1 + ratio * (num_closed_states-1)) / (ratio * mean_ibi),
                "rate close": r_close,
            },
            "activation": {
                "halfmax": 1.,
                "hill exponent": 3.0,
            },
            "ip3": {
                "halfmax": 1,
                "hill exponent": 3.0,
            },
        },
            file, indent=0)


if __name__ == "__main__":
    N = 0
    for i in range(1):
        for j in range(20):
            for k in range(1):
                bool_ca_fix = True
                ca_fix = 0.05*j
                tau = 1 # np.logspace(-1, 2, 50)[j]
                current = 1 # np.logspace(-3, 0, 50)[i]
                bool_adap = False
                tau_adap = 50
                amp_adap = 0.1
                number_cluster =  10
                mean_num_channel = 4
                num_closed_states = 4
                ratio = 0.1
                mean_ibi = 10
                r_close = 50
                ip3 = 1.
                run = k
                write_parameter_file(N, run, bool_ca_fix, ca_fix, tau, current, bool_adap, tau_adap, amp_adap, ip3, number_cluster, mean_num_channel, num_closed_states, ratio, mean_ibi,
                         r_close)
                N+=1
