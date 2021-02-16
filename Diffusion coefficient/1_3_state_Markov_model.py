import numpy as np
from numpy.random import default_rng
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os

def markov_model(ts, dt, r1, r2, r3):
    data = []  # time state
    state = 1
    for t in ts:
        rng = np.random.uniform(0, 1)
        if state == 1:
            pt = r1
        elif state == 2:
            pt = r2
        elif state == 3:
            pt = r3
        if (rng < pt * dt):
            state += 1
            if state == 4:
                state = 1
        data.append([t, state])
    return data

if __name__ == "__main__":
    home = os.path.expanduser("~")
    r1 = 0.1
    r2 = 0.2
    r3 = 0.3

    dt = 0.1
    tmax = 1_000_000
    ts = np.arange(0, tmax, step=dt)
    datas = markov_model(ts, dt, r1, r2, r3)


    file_str = home + "/Data/3_state_Markov_model/3_state_t{:.0e}.dat".format(tmax)
    with open(file_str, "w") as f:
        for data in datas:
            t = data[0]
            x = data[1]
            f.write("{:.1f} {:d} \n".format(t, x))