import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
from matplotlib import rc
from matplotlib import rcParams


def set_default_plot_style():
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = 'Computer Modern Roman'
        rc('text', usetex=True)


def remove_top_right_axis(axis):
        for ax in axis:
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)

set_default_plot_style()
home = os.path.expanduser("~")

CaRest = 0.33
CaFix = 0.4
Ip3 = 1.
n = 4
m = 4

meanI = 10 # for a single channelCaFix
r_opn = 0.13 * np.power(CaFix / CaRest, 3) * (1 + CaRest ** 3) / (1 + CaFix ** 3)
r_ref = 1.3 * np.power(CaFix / CaRest, 3) * (1 + CaRest ** 3) / (1 + CaFix ** 3)
r_cls = 50
pnorm = ((m-1) * r_cls * r_opn + r_cls * r_ref + (n + 1) * (n / 2) * r_opn * r_ref) / (r_cls * r_opn * r_ref)

data = np.loadtxt(home + "/CLionProjects/PhD/calcium_spikes_markov/out/fixed calcium/puff_markov_cafix{:.2f}_tau1.00e+00_j1.00e+00_N10_0.dat".format(CaFix, n))
data_tmp = []
for x in data:
    if x[2] == 0:
        data_tmp.append(x)
data = data_tmp
states = [[], [], [], [], [], [], [], []]
for set1, set2 in zip(data[:-1], data[1:]): #set = [time, state, idx]
    time = set2[0] - set1[0]
    states[int(set1[1])+3].append(time)
probs = [[], [], [], [], [],[], [], []]
for i, s in enumerate(states):
    probs[i] = sum(s)
pmax = sum(probs)
probs = [p/pmax for p in probs]

ptheo =  [1 / r_ref, 1 / r_ref, 1 / r_ref, 1 / r_opn, 1 / r_cls, 2 / r_cls, 3 / r_cls, 4 / r_cls]
ptheo = [p/pnorm for p in ptheo]
for p in ptheo:
    print(p)

popen = ptheo[-4] + ptheo[-3] + ptheo[-2] + ptheo[-1]
def pclosedtheo(i):
    return popen*2*(n+1-i)/(n*(n+1))

plt.plot([1, 2, 3, 4], [pclosedtheo(4), pclosedtheo(3), pclosedtheo(2), pclosedtheo(1)], c="r",zorder=5)
plt.plot([-3, -2, -1, 0, 1, 2, 3, 4], ptheo, c="k", label="Theory")
plt.scatter([-3, -2, -1, 0, 1, 2, 3, 4], [probs[0], probs[1], probs[2] ,probs[3], probs[-1], probs[-2], probs[-3], probs[-4]], label="Simulation \n $\mu=${:.3f}".format(4*probs[-1] + 3*probs[-2] + 2*probs[-3] + probs[-4]))
plt.xticks(ticks=[-3, -2, -1, 0, 1, 2, 3, 4], labels=["$0_4$", "$0_3$", "$0_2$", "$0_1$", "$4$", "$3$", "$2$", "$1$"])
plt.legend()
plt.yscale("log")
plt.savefig(home + "/Data/Calcium/Plots/puff_state_prob_a{:.2f}.pdf".format(CaFix))
plt.show()