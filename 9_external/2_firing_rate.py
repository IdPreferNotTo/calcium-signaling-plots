import subprocess
import numpy as np
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/calcium_firing_rate/cmake-build-release/calcium_firing_rate"

    a = 0.5
    tau = 1.78 #5.62
    j = 0.0355 #0.0126
    ip3s = np.linspace(0.5, 1.5, 101)
    r0s = []
    N = 1

    r0s_file = home + f"/Data/calcium_spikes_theory/ca_langevin_r0_tau{tau:.2e}_j{j:.2e}_over_ip3.dat"
    with open(r0s_file, "w") as f:
        for ip3 in ip3s:
            result = subprocess.run([exetuable, f"{ip3:.2f}", f"{tau:.3f}", f"{j:.3f}", f"{a:.1f}"], stdout=subprocess.PIPE)
            result_string = result.stdout.decode('utf-8')
            print(tau, j, ip3, result_string)
            r0s.append(result_string)
            f.write(f"{ip3:.2e} " + result_string)
            N += 1




