import subprocess
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    exetuable =  home + "/CLionProjects/PhD/calcium/calcium_probability_density/cmake-build-release/calcium_probability_density"
    ip3 = 1.00
    tau = 10.50
    j = 0.011
    a = 0.5
    subprocess.run([exetuable, f"{ip3:.2f}", f"{tau:.3f}", f"{j:.3f}", f"{a:.1f}"])