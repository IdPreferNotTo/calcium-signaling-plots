import numpy as np
from scipy.integrate import quad
import sys
import os

import functions as fc


if __name__ == "__main__":
    #ip3 = int(sys.argv[1])
    #IP3s = np.linspace(0.02, 1, 50)
    IP3s = [1.1]
    for ip3 in IP3s:
        home = os.path.expanduser("~")
        file_str = home + f"/Data/Calcium/data/firing_rate_over_ip3_{ip3:.2f}.dat"
        tau = np.logspace(-1, 2, 50)[33]
        j = np.logspace(-3, 0, 50)[17]

        N = 10
        n = 5
        m = 4
        with open(file_str, "w") as file:
            rate = fc.firing_rate_no_adap(tau, j, N, n, m, ip3)
            file.write(f"{rate:.2e}")
