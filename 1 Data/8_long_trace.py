import os
import numpy as np

if __name__ == "__main__":
    home = os.path.expanduser("~")
    non_stat = [7, 10, 17, 19, 21, 22, 25, 26]
    for j in range(1, 39):
        data = np.loadtxt(home + f"/Data/calcium_spikes_experimental/Spikes/HEK/HEK2/spike_times_{j:d}.dat")
        ieff = data[0]
        T0 = data[1]
        T8 = data[2]
        isis = data[3:]
        isis_without_transient = isis[int(ieff) + 1:]