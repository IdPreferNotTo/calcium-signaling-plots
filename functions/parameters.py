import json
import numpy as np
import os


def write_parameter_file(model, n, i, ip3, bool_ca_fix, ca_fix, tau, p, bool_adap, tau_adap, eps_adap, K, N, M, r_opn_single, r_ref, r_cls, x):
    home = os.path.expanduser('~')
    with open(home + '/CLionProjects/PhD/calcium/calcium_spikes_{}/parameter/{}.json'.format(model, n), 'w') as file:
        json.dump({
            "num parameter": {
                "run": i,
                "output": "BCCN",
                "output puff": False,
                "dt": 10 ** (-4),
                "dt langevin": 10 ** (-3),
                "t_out": 10 ** (-2),
                "max spikes": 5000,
                "max time": 100_000
            },
            "cell": {
                "num cluster": K,
                "timeconstant": tau,
                "blip current": p,
                "ip3": ip3,
                "calcium rest": 0.2,
                "calcium threshold": 0.5,
            },
            "calcium fix":{
                "on":bool_ca_fix,
                "value":ca_fix,
            },
            "adaptation":{
                "on": bool_adap,
                "timeconstant": tau_adap,
                "amplitude": eps_adap,
            },
            "cluster": {
                "num channel": N,
                "number closed states": M,
                "rate open single": r_opn_single,
                "rate ref": r_ref,
                "rate close": r_cls,
            },
            "buffer":{
                "kp": 100,
                "km": 500,
                "bT": x,
            }
        },
            file, indent=0)


if __name__ == "__main__":
    n = 0
    for i in range(61):
        for j in range(61):
            for k in range(1):
                model = "markov"
                ip3 = 1 #np.linspace(0.5, 1.5, 101)[i] #np.linspace(0.02, 2, 100)[i] #IP3 is given realtive the half-max ip3 concentration.
                bool_ca_fix = False
                ca_fix = 0.5
                tau = 14.1 #1.78 #np.logspace(-1, 2, 61)[i]
                p = 0.00794 #0.0355 #np.logspace(-3, 0, 61)[j]
                bool_adap = True
                tau_adap = np.logspace(1, 3, 61)[i]
                eps_adap = np.logspace(-2, 0, 61)[j]
                K = 10
                N = 5
                M = 3
                r_ref = 20
                r_opn_single = 0.1
                r_cls = 50
                run = 5
                x = 0
                write_parameter_file(model, n, run, ip3, bool_ca_fix, ca_fix, tau, p, bool_adap, tau_adap, eps_adap, K, N, M, r_opn_single, r_ref, r_cls, x)
                n+=1
    print(n)
