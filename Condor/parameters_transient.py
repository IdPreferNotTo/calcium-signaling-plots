import json
import numpy as np
import os

def write_parameter_file(n, i, ip3, bool_ca_fix, ca_fix, tau, p, bool_adap, tau_adap, eps_adap, K, N, M, r_opn_single, r_ref, r_cls, x):
    home = os.path.expanduser('~')
    with open(home + f'/CLionProjects/PhD/calcium/calcium_spikes_markov_transient/parameter/{n:d}.json', 'w') as file:
        json.dump({
            "num parameter": {
                "run": i,
                "output": "local",
                "output puff": True,
                "dt": 10 ** (-3),
                "dt langevin": 10 ** (-2),
                "t_out": 10 ** (-2),
                "max spikes": 5000,
                "max time": 1_000_000
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
                "kp": 10,
                "km": 50,
                "bT": x,
            }
        },
            file, indent=0)


if __name__ == "__main__":
    n0 = 63
    n = 0
    for i in range(21):
        for j in range(1):
            for k in range(1):
                K = 10
                N = 5
                M = 3
                r_ref = 20
                r_opn_single = 0.1
                r_cls = 50
                ip3 = 1
                bool_ca_fix = False
                ca_fix = 0.2
                tau = 1
                p = 0.06
                bool_adap = True
                tau_adap = 200 # np.logspace(1, 3, 21)[i]
                eps_adap = np.logspace(-2, 0, 21)[i] # 0.05
                print(i+n0, f"{tau_adap:.2e}")
                run = 5
                x = 0
                write_parameter_file(n0 + n, run, ip3, bool_ca_fix, ca_fix, tau, p, bool_adap, tau_adap, eps_adap, K, N, M, r_opn_single, r_ref, r_cls, x)
                n+=1
    print(n0 + n)
