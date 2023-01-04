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
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin/cmake-build-release/calcium_spikes_langevin"

    a = 0.5
    model = "langevin"
    ip3 = 1
    tau = 1
    p = 0.064

    r_cls = 50
    r_ref = 20
    r_opn_single = 0.1
    ropnmax = 5 * r_opn_single * ((1. + np.power(0.20, 3)) / np.power(0.20, 3))
    r_opn_ct = ropnmax * (np.power(0.5, 3) / (1. + np.power(0.5, 3)))
    mean_puff = (6) * (7) / (6 * r_cls)
    tau_tot = 1 / r_opn_ct + 2 / r_ref + 6 / (2 * r_cls)
    mean_puff_ct = mean_puff / tau_tot
    p_bif  = ((0.5 - 0.2) / (10 * mean_puff_ct * tau))

    print(p_bif)
    K = 10
    N = 5
    M = 3
    r_ref = 20
    r_opn_single = 0.1
    r_cls = 50

    write_parameter_file(model, 0, 0, ip3, False, 0.5, tau, p, False, 0, 0, K, N, M, r_opn_single, r_ref, r_cls, 0)
    subprocess.run([exetuable, "0"], stdout=subprocess.PIPE)

    home = os.path.expanduser("~")
    folder_langevin = home + "/CLionProjects/PhD/calcium/calcium_spikes_langevin/out/"
    file_spikes_langevin = folder_langevin + f"spike_times_langevin_ip1.00_tau{tau:.2e}_j{p:.2e}_K10_0.dat"
    ISIs = np.loadtxt(file_spikes_langevin)
    mean_ISIs = np.mean(ISIs)
    cv_ISIs = np.std(ISIs)/mean_ISIs
    print("mean: ", mean_ISIs, "cv: ", cv_ISIs)






