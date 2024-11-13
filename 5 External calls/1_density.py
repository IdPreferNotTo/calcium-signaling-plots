import subprocess
import os

if __name__ == "__main__":
    home = os.path.expanduser("~")
    exetuable =  home + "/CLionProjects/PhD/calcium/theory_probability_density/cmake-build-release/theory_probability_density"

    tau = 5
    j = 0.015
    cer = 1.0
    subprocess.run([exetuable, f"{tau:.3f}", f"{j:.3f}", f"{cer:.1f}"])