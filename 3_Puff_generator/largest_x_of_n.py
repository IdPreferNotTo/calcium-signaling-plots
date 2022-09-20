import numpy as np
import random
import matplotlib.pyplot as plt

def fac(n):
    fac = 1
    for i in range(1, n+1):
        fac *= i
    return fac

if __name__ == "__main__":
    rate = 0.3
    n = 10
    max_Xs = []
    for i in range(100_000):
        Xs = np.random.exponential(1./rate, n)
        max_Xs.append(max(Xs))

    ts = np.linspace(0, 50, 500)
    ys = rate*n*(1-np.exp(-rate*ts))**(n-1)*np.exp(-rate*ts)

    rate_eff = 4*rate
    ys2 = np.exp(-rate_eff*ts)*(rate_eff**n)*ts**(n-1)/fac(n-1)

    plt.plot(ts, ys, c="k")
    plt.plot(ts, ys2, c="k", ls=":")
    plt.hist(max_Xs, bins=100, density=True)
    plt.show()

