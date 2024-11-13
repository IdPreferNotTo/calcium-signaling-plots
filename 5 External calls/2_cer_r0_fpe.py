import subprocess
import numpy as np
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    exetuable = home + "/CLionProjects/PhD/calcium/theory_firing_rate/cmake-build-release/theory_firing_rate"

    tau = 5
    j = 0.015
    r0s_file = home + f"/Data/calcium/theory/cer_r0_fpe_tau{tau:.2e}_j{j:.2e}.dat"
    cers = np.linspace(0.5, 1.0, 401)
    with open(r0s_file, "w") as f:
        for cer in cers:
            print(cer)
            result = subprocess.run([exetuable, f"{tau:.8f}", f"{j:.8f}", f"{cer:.8f}"], stdout=subprocess.PIPE)
            result_string = result.stdout.decode('utf-8')
            r0_self = float(result_string)
            f.write(f"{cer:.4e} {r0_self:.4e} \n")


    tau = 1
    j = 0.06
    r0s_file = home + f"/Data/calcium/theory/cer_r0_fpe_tau{tau:.2e}_j{j:.2e}.dat"
    cers = np.linspace(0.5, 1.0, 401)
    with open(r0s_file, "w") as f:
        for cer in cers:
            print(cer)
            result = subprocess.run([exetuable, f"{tau:.8f}", f"{j:.8f}", f"{cer:.8f}"], stdout=subprocess.PIPE)
            result_string = result.stdout.decode('utf-8')
            r0_self = float(result_string)
            f.write(f"{cer:.4e} {r0_self:.4e} \n")








