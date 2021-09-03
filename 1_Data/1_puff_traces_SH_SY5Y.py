import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

home = os.path.expanduser("~")
file_str = home + "/Desktop/Ca data/Puffs/SH_SY5Y_puffs_longtraces.dat"
data = np.loadtxt(file_str)

n = len(data[0])
for j in range(1, n):
    row = [x[j] for x in data]
    idxs = [i for (i, x) in enumerate(row) if x > 500]
    print(idxs)
    for i, idx in enumerate(idxs):
        data = np.delete(data, (idx-i), axis=0)

    n_avr = 50
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
    gs = gridspec.GridSpec(3, 1)
    ax2 = fig.add_subplot(gs[1:3, 0])
    ax1 = fig.add_subplot(gs[0, 0], sharex=ax2)
    axis = [ax1, ax2]
    axin = inset_axes(ax2, width="40%", height="40%", loc=1)
    ipis = []
    for t1, t2 in zip(puff_times[:-1], puff_times[1:]):
        ipis.append(t2-t1)
    ipi_mean = np.mean(ipis)
    ipi_var =  np.var(ipis)
    Cv = np.sqrt(ipi_var)/ipi_mean

    #cas_unbias = [ca - avr_ca for ca, avr_ca in zip(cas, avr_cas)]
    for ptime in puff_times:
        ax1.axvline(ptime)
    ax2.plot(times[n_avr:-n_avr], cas[n_avr:-n_avr])
    ax2.plot(times[n_avr:-n_avr], avr_cas, c="k")
    axin.hist(ipis, density=True)
    inv_gauss = []

    ts = np.linspace(0, 20, 200)
    for t in ts[1:]:
        pt = (ipi_mean/t)**(3/2)*(1/(np.sqrt(2*np.pi)*Cv*ipi_mean))*np.exp(-(ipi_mean*(t -ipi_mean)**2)/(2*t*(Cv*ipi_mean)**2))
        inv_gauss.append(pt)
    axin.plot(ts[1:], inv_gauss, c="C3")
    #ax2.plot(times[n_avr:-n_avr], [x + 20 for x in avr_cas], c="C7")
    #ax2.set_ylim([80, 200])
    #ax2.set_xlim([0,50])
    plt.savefig(home + "/Desktop/Ca data/Puffs/Plots/SH_{:d}.png".format(j), transparent=True)
    plt.show()
    plt.close()