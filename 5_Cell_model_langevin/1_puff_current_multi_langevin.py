import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

def deterministic_integrate_multi_langein(n01, n02, n2, n1, rref, ropn, rcls, dt):
    dn01 = (n1*rcls - n01*rref)*dt
    dn02 = (n01*rref - n02*ropn)*dt
    dn2 = (0.5*n02*ropn - n2*rcls)*dt
    dn1 = (0.5*n02*ropn + n2*rcls - n1*rcls)*dt
    n01 += dn01
    n02 += dn02
    n2 += dn2
    n1 += dn1
    return n01, n02, n2, n1


def integrate_multi_langein(n01, n02, n2, n1, rref, ropn, rcls, dt):
    xi01 = np.random.normal(0, 1)
    xi02 = np.random.normal(0, 1)
    xi2 = np.random.normal(0, 1)
    xi1 = np.random.normal(0, 1)
    noise01 = np.sqrt(abs(n01)*rref*dt)*xi01
    noise02 = np.sqrt(abs(n02)*ropn*dt)*xi02
    noise2 = np.sqrt(abs(n2)*rcls*dt)*xi2
    noise1 = np.sqrt(abs(n1)*rcls*dt)*xi1
    dn01 = (n1*rcls - n01*rref)*dt - noise1 + noise01
    dn02 = (n01*rref - n02*ropn)*dt - noise01 + noise02 #+ np.sqrt(0.5*abs(n02)*ropn*dt)*xi022
    dn2 = (0.5*n02*ropn - n2*rcls)*dt - noise02/2 + noise2
    dn1 = (0.5*n02*ropn + n2*rcls - n1*rcls)*dt - noise02/2 - noise2 + noise1
    n01 += dn01
    n02 += dn02
    n2 += dn2
    n1 += dn1
    return n01, n02, n2, n1

if __name__ == "__main__":
    start = 0
    stop = 2
    steps = 100000
    dt = (stop-start)/steps
    ts =  np.linspace(start, stop, steps)
    rref = 5.2
    ropn = 0.52
    rcls = 50
    N = 10
    n01 = 2.5
    n02 = 2.5
    n2 = 2.5
    n1 = 2.5
    nss = []
    ns = [n01, n02, n2, n1]
    for t in ts[:int(steps/2)]:
        ns = deterministic_integrate_multi_langein(*ns, rref, ropn, rcls, dt)
        nss.append(ns)
    for t in ts[int(steps/2):]:
        ns = integrate_multi_langein(*ns, rref, ropn, rcls, dt)
        nss.append(ns)

    n01s, n02s, n2s, n1s = np.transpose(nss)
    fig = plt.figure(tight_layout=True, figsize=(4, 5))
    subplots = 6
    gs = gs.GridSpec(subplots, 1)
    axis = []
    for i in range(subplots):
        ax = fig.add_subplot(gs[i, 0])
        axis.append(ax)
    for i in range(subplots):
        if i==(subplots-2):
            axis[i].plot(ts, [x[0] + x[1] for x in nss], c="C1")
            axis[i].set_xticks([])
        elif i==(subplots-1):
            axis[i].plot(ts, [x[2] + x[3] for x in nss], c="C1")
        else:
            axis[i].plot(ts, [x[i] for x in nss])
            axis[i].set_xticks([])
    plt.show()
