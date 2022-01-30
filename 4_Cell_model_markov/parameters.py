import json
import numpy as np
import os


def write_parameter_file(model, N, i, ip3, bool_ca_fix, ca_fix, tau, current, bool_adap, tau_adap, amp_adap, number_cluster, num_channel, num_closed_states, ratio, mean_ibi, r_close):
    home = os.path.expanduser('~')
    with open(home + '/CLionProjects/PhD/calcium_spikes_{}/parameter/{}.json'.format(model, N), 'w') as file:
        json.dump({
            "num parameter": {
                "run": i,
                "output": "BCCN",
                "dt": 10 ** (-4),
                "dt langevin": 10 ** (-2),
                "t_out": 10 ** (-1),
                "max spikes": 10_000,
                "max time": 1_000_000
            },
            "cell": {
                "num cluster": number_cluster,
                "timeconstant": tau,
                "blip current": current,
                "ip3": ip3,
                "calcium rest": 0.33,
                "calcium threshold": 1.0,
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
                "num channel": num_channel,
                "number closed states": num_closed_states,
                "rate open single": (1. + ratio * (num_closed_states-1.)) / mean_ibi,
                "rate ref single": (1. + ratio * (num_closed_states-1.)) / (ratio * mean_ibi),
                "rate close": r_close,
            },
        },
            file, indent=0)


if __name__ == "__main__":
    N = 0
    for i in range(100):
        for j in range(1):
            for k in range(1):
                print(N)
                model = "markov"
                ip3 = 1 #IP3 is given realtive the half-max ip3 concentration.
                bool_ca_fix = False
                ca_fix = 0.33
                tau = 2.0 #np.logspace(-1, 2, 31)[i]
                current = 0.2*((i+1)/100) #np.logspace(-3, 0, 31)[j]
                bool_adap = True
                tau_adap = 100 #np.logspace(1, 3, 41)[i] #300
                amp_adap = 0.05 #np.logspace(-2, 0, 41)[j] #0.1
                number_cluster = 10
                num_channel = 5
                num_closed_states = 4
                ratio = 0.1
                mean_ibi = 10
                r_close = 50
                run = 0
                write_parameter_file(model, N, run, ip3, bool_ca_fix, ca_fix, tau, current, bool_adap, tau_adap, amp_adap, number_cluster, num_channel, num_closed_states, ratio, mean_ibi,
                         r_close)
                N+=1
