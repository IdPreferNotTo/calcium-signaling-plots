import subprocess
import numpy as np
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/calcium_firing_rate/cmake-build-release/calcium_firing_rate"

    parameter_set = 1
    #taus = [5.0, 1.0]
    #ps = [0.015, 0.06]
    #eps_er_fix = 0.1
    #tau_er_fix = 100
    #tau_ers = np.logspace(1, 3, 21)
    #eps_ers = np.logspace(-2, 0, 21)

    if parameter_set == 1:
        tau = 5
        j = 0.015
        tau_er = 100
        eps_ers = np.logspace(-2, 0, 21)
        r0s = np.linspace(0, 0.05, 101)

    for eps_er in eps_ers:
        r0s_file = home + f"/Data/calcium_spikes_theory/r0_selfcon_tau{tau:.2e}_j{j:.2e}_{eps_er:.2e}_{tau_er:.2e}.dat"
        with open(r0s_file, "w") as f:
            for r0 in r0s:
                print(r0)
                result = subprocess.run([exetuable, f"{1.:.2f}", f"{tau:.3f}", f"{j:.3f}", f"{eps_er:.2f}" , f"{tau_er:.1f}", f"{r0:.8f}"], stdout=subprocess.PIPE)
                result_string = result.stdout.decode('utf-8')
                r0_self = float(result_string)
                f.write(f"{r0:.2e} {r0_self:.2e} \n")




