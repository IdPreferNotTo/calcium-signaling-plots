import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut

def post_puff_ca(ts, cas, starts, stops):
    if stops[-1] < starts[-1]:
        del starts[-1]
    i = 1
    t_post_puff_averages = []
    ca_post_puff_averages = []
    t_post_puff_average = []
    ca_post_puff_average = []
    for t, ca in zip(reversed(ts), reversed(cas)):
        if t < stops[-i] and t > starts[-i]:
            ca_post_puff_average.append(ca)
            t_post_puff_average.append(t - starts[-i])
        if t < starts[-i]:
            ca_post_puff_averages.append(list(ca_post_puff_average))
            t_post_puff_averages.append(list(t_post_puff_average))
            ca_post_puff_average.clear()
            t_post_puff_average.clear()
            i+=1
        if i == len(starts):
            break
    return [t_post_puff_averages, ca_post_puff_averages]


def pre_puff_ca(ts, cas, starts, stops):
    if stops[-1] > starts[-1]:
        del stops[-1]
    i = 1
    t_puff_triggered_averages = []
    ca_puff_triggered_averages = []
    t_puff_triggered_average = []
    ca_puff_triggered_average = []
    for t, ca in zip(reversed(ts), reversed(cas)):
        if t > stops[-i] and t < starts[-i]:
            ca_puff_triggered_average.append(ca)
            t_puff_triggered_average.append(t - starts[-i])
        if t < stops[-i]:
            ca_puff_triggered_averages.append(list(ca_puff_triggered_average))
            t_puff_triggered_averages.append(list(t_puff_triggered_average))
            ca_puff_triggered_average.clear()
            t_puff_triggered_average.clear()
            i+=1
        if i == len(starts):
            break
    return [t_puff_triggered_averages, ca_puff_triggered_averages]


def pre_puff_average(tss, cass):
    n_max = max([len(cas) for cas in cass])
    n_most = int(0.7*n_max)
    ca_mean_pre = []
    t_mean_pre = []
    for n in range(n_most):
        cas_trans= []
        for ts, cas in zip(tss, cass):
            if len(cas) > n:
                cas_trans.append(cas[n])
                t = ts[n]
        ca_mean_pre.append(np.mean(cas_trans))
        t_mean_pre.append(t)
    return t_mean_pre, ca_mean_pre


def post_puff_average(tss, cass):
    n_max = max([len(cas) for cas in cass])
    n_most = int(0.5*n_max)
    ca_mean_post = []
    t_mean_post = []
    for n in range(1, n_most):
        cas_trans= []
        for ts, cas in zip(tss, cass):
            if len(cas) > n:
                cas_trans.append(cas[-n])
                t = ts[-n]
        ca_mean_post.append(np.mean(cas_trans))
        t_mean_post.append(t)
    return t_mean_post, ca_mean_post


def derivation_ca(ts, cas):
    x = []
    dcas = []
    dt = ts[1]-ts[0]
    for ca0, ca1 in zip(cas[:-1], cas[1:]):
        #remember that cas is backwards thus ca1-ca2 instead of ca2-ca1
        dcas.append((ca1 - ca0)/dt)
        x.append(ca0)
    return x, dcas


f, ax = plt.subplots(1, 1, figsize=(6, 9 / 2))
bins = 100
ip3 = 0.3
N = 64
data = np.loadtxt("/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(N, ip3), skiprows=10_000)
ts, cas, j1, n_open = np.transpose(data)
ipis, starts, stops = ut.spike_times(ts, cas, n_open, ip3)
tca_puff_triggered_averages = pre_puff_ca(ts, cas, starts, stops)
tca_puff_followed_averages = post_puff_ca(ts, cas, starts, stops)

ts_pre, cas_pre = tca_puff_triggered_averages
ts_post, cas_post = tca_puff_followed_averages

t_mean_pre, ca_mean_pre = pre_puff_average(ts_pre, cas_pre)
t_mean_post, ca_mean_post = post_puff_average(ts_post, cas_post)
xs, dca_mean_pre = derivation_ca(t_mean_pre, ca_mean_pre)
#ax.plot(xs, dca_mean_pre, c="k")

minimum = 1
maximum = 0
for ca_pre in cas_pre:
    if len(ca_pre) > 0:
        if maximum < max(ca_pre):
            maximum = max(ca_pre)
        if minimum > min(ca_pre):
            minimum = min(ca_pre)


xs_pre_averaged = [minimum + i*(maximum-minimum)/bins for i in range(bins)]
dcas_pre_all= empty_lists = [[] for _ in range(bins)]


for n, (ts, cas) in enumerate(zip(ts_pre, cas_pre)):
    if len(cas) > 2:
        xs, dcas_pre = derivation_ca(ts, cas)
        if n < 100:
            plt.scatter(xs, dcas_pre, c="C0", alpha=0.1, s=1)
        for x, dca in zip(xs, dcas_pre):
            index = int((x - minimum)*bins/(maximum - minimum))
            dcas_pre_all[index-1].append(dca)

dcas_pre_averaged = []
for dcas_pre in dcas_pre_all:
    dcas_pre_averaged.append(np.mean(dcas_pre))

plt.plot(xs_pre_averaged, dcas_pre_averaged, c="k")
plt.hlines(0, minimum, maximum, ls="--", zorder=3)

fit1 = []
fit2 = []
fit3 = []
ca_rest = 0.09
ca_th = 0.08
tau = 20
Delta = 0.02
cas = np.linspace(minimum, maximum, 200)
for ca in cas:
    fit1.append(-(ca - ca_rest)/tau)
    fit2.append((Delta/tau) *np.exp((ca - ca_th)/Delta))
    fit3.append(-(ca - ca_rest)/tau + (Delta/tau) *np.exp((ca - ca_th)/Delta))
plt.plot(cas, fit1, c="C3", ls = "--")
plt.plot(cas, fit2, c="C3", ls = "--")
plt.plot(cas, fit3, c="C3", label=r"Exp. IF: $\tau \dot{x}= -(x-x_R) + \Delta$exp$[(x - x_T)/\Delta]$" +"\n" r"$\tau = {} x_R = {}, x_T = {}, \Delta = {},$".format(tau, ca_rest, ca_th, Delta))
plt.legend()
ax.set_ylim([-0.01, 0.05])
ax.set_ylabel(r"$\Delta$Ca$^{2+}/dt$")
ax.set_xlabel(r"Ca$^{2+}$")
plt.savefig("/home/lukas/Plots/lr_ca_phd/lr_pre_spike_dca_ip{:.2f}_N{:d}.pdf".format(ip3, N), transparent=True)
plt.show()