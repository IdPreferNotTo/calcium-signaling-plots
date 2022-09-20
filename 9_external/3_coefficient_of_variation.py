import subprocess
import numpy as np
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/calcium_coefficient_of_variation/cmake-build-release/calcium_coefficient_of_variation"
    ip3 = 1.00
    a = 0.5
    if a == 0.0:
        interpretation = "ito"
    elif a == 0.5:
        interpretation = "strat"
    else:
        interpretation = "haenngi"
    tau = 4.47  # [4.47, 1.26]
    j = 0.0178  # [0.0178, 0.0562]
    ip3s = np.linspace(0.02, 2, 100)
    r0s = []
    N = 1

    for ip3 in ip3s:
        result = subprocess.run([exetuable, f"{ip3:.2f}", f"{tau:.3f}", f"{j:.3f}", f"{a:.1f}"],
                                stdout=subprocess.PIPE)
        result_string = result.stdout.decode('utf-8')
        print(N, result_string)
        r0s.append(result_string)
        N += 1

    r0s_file = home + f"/Data/calcium_spikes_theory/ca_langevin_{interpretation:s}_cvs.dat"
    with open(r0s_file, "w") as f:
        for ip3, r0 in zip(ip3s, r0s):
            f.write(f"{ip3:.2f} " + r0)
