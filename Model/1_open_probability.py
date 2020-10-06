import numpy as np
import matplotlib.pyplot as plt

p_max = 0.8
Kact = 0.21
Hact = 1.9
Hinh = 3.9
Kinf = 52
Kip3 = 0.05
Hip3 = 4
ip3 = 0.1
Kinh = Kinf*np.power(1+np.power(Kip3/ip3, Hinh), -1)

xs = np.logspace(-2, 2, 100)
ps1 = []
ps2 = []
ps3 = []
for x in xs:
    ps1.append(p_max*np.power(1+np.power(Kact/x, Hact), -1)*np.power(1 + np.power(x/Kinh, Hinh), -1))
    ps2.append(p_max*np.power(1+np.power(Kact/x, Hact), -1))
    ps3.append(p_max*np.power(1 + np.power(x/Kinh, Hinh), -1))
plt.plot(xs, ps1, c="k")
plt.plot(xs, ps2, alpha=0.7)
plt.plot(xs, ps3, alpha=0.7)
plt.xscale("log")
plt.show()
