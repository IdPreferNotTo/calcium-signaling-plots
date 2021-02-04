import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

home = os.path.expanduser("~")

CaRest = 0.33
CaFix = 0.35
Ip3 = 1.
NCha = 4
NClosedStates = 4
m = NClosedStates - 2
meanI = 10 # for a single channelCaFix
ropen = 0.13*np.power(CaFix/CaRest, 3)*(1+CaRest**3)/(1+CaFix**3)
rref = 1.3*np.power(CaFix/CaRest, 3)*(1+CaRest**3)/(1+CaFix**3)
rclose = 50
pnorm = (2*(m+1)*rclose*ropen + 2*rclose*rref + (NCha+1)*NCha*ropen*rref)/(2*NCha*ropen*rref)

data = np.loadtxt(home + "/CLionProjects/calcium-spikes-from-puff-phd/out/fixed calcium/puff_fixedCa{:.2f}_tau1.00e+00_j1.00e+00_Ncha{:d}_Nclu10_Ncls4_rO0.13_rC50.00_rR1.30_N0.dat".format(CaFix, NCha))
data_tmp = []
for x in data:
    if x[2] == 3:
        data_tmp.append(x)
data = data_tmp
states = [[], [], [], [], [],[], [], []]
for set1, set2 in zip(data[:-1], data[1:]): #set = [time, state, idx]
    time = set2[0] - set1[0]
    states[int(set1[1])+3].append(time)
probs = [[], [], [], [], [],[], [], []]
for i, s in enumerate(states):
    probs[i] = sum(s)
for p in probs:
    print(p)

pmax = sum(probs)
probs = [p/pmax for p in probs]
print(sum(probs))
ptheo =  [rclose/(NCha*rref), rclose/(NCha*rref), rclose/(NCha*rref), rclose/(NCha*ropen), 1/NCha, 2/NCha, 3/NCha, 4/NCha]
ptheo = [p/pnorm for p in ptheo]

popen = ptheo[-4] + ptheo[-3] + ptheo[-2] + ptheo[-1]
def pclosedtheo(i):
    return popen*2*(NCha+1-i)/(NCha*(NCha+1))

plt.plot([1, 2, 3, 4], [pclosedtheo(4), pclosedtheo(3), pclosedtheo(2), pclosedtheo(1)], c="r",zorder=5)
plt.plot([-3, -2, -1, 0, 1, 2, 3, 4], ptheo, c="k", label="Theory")
plt.scatter([-3, -2, -1, 0, 1, 2, 3, 4], [probs[0], probs[1], probs[2] ,probs[3], probs[-1], probs[-2], probs[-3], probs[-4]], label="Simulation")
plt.xticks(ticks=[-3, -2, -1, 0, 1, 2, 3, 4], labels=["$0_4$", "$0_3$", "$0_2$", "$0_1$", "$4$", "$3$", "$2$", "$1$"])
plt.legend()
plt.yscale("log")
plt.savefig(home + "/Data/Calcium/Plots/puff_state_prob_fixedCa{:.2f}.pdf".format(CaFix))
plt.show()