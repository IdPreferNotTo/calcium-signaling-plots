import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

home = os.path.expanduser("~")
#file_str = home + "/Desktop/puff_data/puffs_thurley_t1.dat"
file_str = home + "/Desktop/Ca data/SH_SY5Y_puffs_longtraces.dat"
#file_str = home + "/Desktop/puff_data/HEK_traces/0712runb6_2"
data = np.loadtxt(file_str)

n = len(data[0])
for j in range(1, n):
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    print(idxs)
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx-i), axis=0)

    n_avr = 100
    stat_cas = []
    times = [x[0]/1000 for x in data]
    cas = [x[j] for x in data]
    puff_times =[]
    avr_cas = []
    puff = False

    for idx, (t, y) in enumerate(zip(times, cas)):
        if idx < n_avr or idx >= (len(cas)-n_avr):
            continue
        mov_avr = np.mean(cas[idx-n_avr: idx])
        mov_std = np.std(cas[idx-n_avr: idx])
        avr_cas.append(mov_avr)
        stat_cas.append(y-mov_avr)
        if y > 1.2*mov_avr and puff == False:
            puff_times.append(t)
            puff = True
        if y < mov_avr + 20 and puff == True:
            puff = False

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    ipis = []
    for t1, t2 in zip(puff_times[:-1], puff_times[1:]):
        ipis.append(t2-t1)
    ipi_mean = np.mean(ipis)
    ipi_var =  np.var(ipis)
    Cv = np.sqrt(ipi_var)/ipi_mean

    ts = np.linspace(0, 20, 200)
    inv_gauss = []
    for t in ts[1:]:
        pt = (ipi_mean/t)**(3/2)*(1/(np.sqrt(2*np.pi)*Cv*ipi_mean))*np.exp(-(ipi_mean*(t -ipi_mean)**2)/(2*t*(Cv*ipi_mean)**2))
        inv_gauss.append(pt)
    ax.plot(ts[1:], inv_gauss, c="C3")
    ax.hist(ipis, density=True)
    plt.show()
    plt.close()
