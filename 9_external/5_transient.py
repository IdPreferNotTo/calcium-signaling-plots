import subprocess
import numpy as np
import json
import os

def write_parameter_file(model, n, i, ip3, bool_ca_fix, ca_fix, tau, p, bool_adap, tau_adap, eps_adap, K, N, M, r_opn_single, r_ref, r_cls, x):
    home = os.path.expanduser('~')
    with open(home + '/CLionProjects/PhD/calcium/calcium_spikes_{}/parameter/{}.json'.format(model, n), 'w') as file:
        json.dump({
            "num parameter": {
                "run": i,
                "output": "local",
                "output puff": False,
                "dt": 10 ** (-4),
                "dt langevin": 10 ** (-3),
                "t_out": 10 ** (-2),
                "max spikes": 25,
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
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin/cmake-build-release/calcium_spikes_langevin"

    a = 0.5
    model = "markov"
    ip3 = 1
    tau = 5
    p = 0.015
    tau_ers = np.logspace(1, 3, 41)
    eps_ers = np.logspace(-2, 0, 41)
    K = 10
    N = 5
    M = 3
    r_ref = 20
    r_opn_single = 0.1
    r_cls = 50
    for tau_er in tau_ers[0:1]:
        for eps_er in eps_ers[0:1]:
            ISIss = []
            for i in range(100):
                print(i)
                write_parameter_file(model, 0, 0, ip3, False, 0.5, tau, p, True, tau_er, eps_er, K, N, M, r_opn_single, r_ref, r_cls, 0)
                subprocess.run([exetuable, "0"], stdout=subprocess.PIPE)

                home = os.path.expanduser("~")
                folder_langevin = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin/out/"
                file_spikes_langevin = folder_langevin + f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{p:.2e}_K10_0.dat"
                ISIs = np.loadtxt(file_spikes_langevin)
                ISIss.append(ISIs)







