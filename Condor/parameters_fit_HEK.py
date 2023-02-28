import json
import numpy as np
import os

def write_parameter_file(model, n, i, ip3, bool_ca_fix, ca_fix, tau, p, bool_adap, tau_adap, eps_adap, K, N, M, r_opn_single, r_ref, r_cls, x):
    home = os.path.expanduser('~')
    with open(home + '/CLionProjects/PhD/calcium/calcium_spikes_{}/parameter/{}.json'.format(model, n), 'w') as file:
        json.dump({
            "num parameter": {
                "run": i,
                "output": "local",
                "output puff": True,
                "dt": 10 ** (-3),
                "dt langevin": 10 ** (-2),
                "t_out": 10 ** (-2),
                "max spikes": 5_000,
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
    taus = [10.28, 13.38, 11.37,
            4.532, 4.180, 9.692,
            2.741, 5.492, 5.492,
            5.492, 65., 6.537,
            11.27, 2.407, 12.08,
            2.781, 3.988, 10.27,
            18.74, 18.74, 7.619,
            15.49, 7.777, 6.850]

    ps = [8.779e-03, 9.923e-03, 1.168e-02,
          2.288e-02, 1.866e-02, 1.308e-02,
          2.510e-02, 1.471e-02, 1.471e-02,
          1.471e-02, 6.70e-03, 1.501e-02,
          1.139e-02, 3.133e-01, 8.519e-03,
          2.666e-02, 2.029e-02, 8.700e-03,
          6.000e-03, 6.000e-03, 1.163e-02,
          7.463e-03, 1.087e-02, 1.749e-02]

    tau_ers = [4.780e+02, 7.934e+02, 2.365e+03,
               6.114e+02, 7.051e+02, 1.026e+03,
               1.843e+02, 8.263e+02, 5.111e+02,
               1.046e+03, 1.071e+03, 6.345e+02,
               1.123e+03, 2.197e+02, 1.761e+02,
               3.066e+02, 1.668e+03, 4.249e+02,
               1.067e+03, 1.052e+03, 3.837e+02,
               2.251e+03, 5.702e+02, 1.395e+03]

    eps_ers = [3.472e-01, 1.294e-01, 9.794e-02,
               2.512e-01, 6.998e-02, 2.801e-01,
               1.628e-01, 1.202e-01, 1.918e-01,
               3.848e-02, 2.965e-01, 1.128e-01,
               1.318e-01, 1.356e-01, 2.065e-01,
               1.237e-01, 8.027e-02, 1.162e-01,
               1.108e-01, 1.130e-01, 1.797e-01,
               7.569e-02, 5.251e-02, 1.093e-01]
    n = 0
    for i in range(24):
        for j in range(1):
            for k in range(1):
                model = "markov"
                K = 10
                N = 5
                M = 3
                r_ref = 20
                r_opn_single = 0.1
                r_cls = 50
                ip3 = 1
                bool_ca_fix = False
                ca_fix = 0.2
                tau = taus[i]
                p = ps[i]
                bool_adap = True
                tau_adap = tau_ers[i]
                eps_adap = eps_ers[i]
                run = 5
                x = 0
                write_parameter_file(model, n, run, ip3, bool_ca_fix, ca_fix, tau, p, bool_adap, tau_adap, eps_adap, K, N, M, r_opn_single, r_ref, r_cls, x)
                n+=1
    print(n)
