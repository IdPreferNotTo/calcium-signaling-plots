import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def post_puff_ca(ts, cas, starts, stops):
    if stops[-1] < starts[-1]:
        del starts[-1]
    i = 1
    t_puff_following_averages = []
    ca_puff_following_averages = []
    t_puff_following_average = []
    ca_puff_following_average = []
    for t, ca in zip(reversed(ts), reversed(cas)):
        if t < stops[-i] and t > starts[-i]:
            ca_puff_following_average.append(ca)
            t_puff_following_average.append(t - starts[-i])
        if t < starts[-i]:
            ca_puff_following_averages.append(list(ca_puff_following_average))
            t_puff_following_averages.append(list(t_puff_following_average))
            ca_puff_following_average.clear()
            t_puff_following_average.clear()
            i+=1
        if i == len(starts):
            break
    return [t_puff_following_averages, ca_puff_following_averages]


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
    ca_mean_pre = []
    t_mean_pre = []
    for n in range(n_max):
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
    ca_mean_post = []
    t_mean_post = []
    for n in range(1, n_max):
        cas_trans= []
        for ts, cas in zip(tss, cass):
            if len(cas) > n:
                cas_trans.append(cas[-n])
                t = ts[-n]
        ca_mean_post.append(np.mean(cas_trans))
        t_mean_post.append(t)
    return t_mean_post, ca_mean_post


def derivation_ca(ts, cas):
    xs = []
    dcas = []
    dt = ts[1]-ts[0]
    for t, ca0, ca1 in zip(ts, cas[:-1], cas[1:]):
        if t > -10 and t < 2:
            dcas.append((ca1 - ca0)/dt)
            xs.append(ca0)
    return xs, dcas


f, ax = plt.subplots(1, 1, figsize=(6, 9 / 2))
ip3 = 0.2
N = 16
data = np.loadtxt("/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(N, ip3), skiprows=10_000)
ts, cas, j1, n_open = np.transpose(data)
ipis, starts, stops = ut.spike_times(ts, cas, n_open, ip3)
tca_puff_triggered_averages = pre_puff_ca(ts, cas, starts, stops)
tca_puff_followed_averages = post_puff_ca(ts, cas, starts, stops)

ts_pre, cas_pre = tca_puff_triggered_averages
ts_post, cas_post = tca_puff_followed_averages
for t_pre, ca_pre, t_post, ca_post in zip(ts_pre, cas_pre, ts_post, cas_post):
    ax.plot(t_pre, ca_pre, c="C0", alpha=0.2)
    ax.plot(t_post, ca_post, c="C1", alpha=0.2)

t_mean_pre, ca_mean_pre = pre_puff_average(ts_pre, cas_pre)
t_mean_post, ca_mean_post = post_puff_average(ts_post, cas_post)

ax.plot(t_mean_pre, ca_mean_pre, c="k", label="$N_{{ch}} = {:d}$".format(N))
ax.plot(t_mean_post, ca_mean_post, c="k", label="$N_{{ch}} = {:d}$".format(N))
ax.legend()

xs1, dcas1 = derivation_ca(t_mean_pre, ca_mean_pre)
xs2, dcas2 = derivation_ca(t_mean_post, ca_mean_post)
ax.axvline(-10, ls="--", c="C7")
ax.axvline(2, ls="--", c="C7")
ax_inset = inset_axes(ax, width=3, height=1, loc=1)
ax_inset.plot(xs1, dcas1)
ax_inset.plot(xs2, dcas2)
ax_inset.set_xlabel("Ca$^{2+}$")
ax_inset.set_ylabel(r"$\Delta$Ca$^{2+}/dt$")
ca0, h0 = ut.fixed_point(0.3)
print(ca0, h0)
ax.set_xlim([-20, 10])
ax.set_ylabel("Ca$^{2+}$")
ax.set_xlabel(r"$t_{rel}$")
plt.savefig("/home/lukas/Plots/lr_ca_phd/lr_pre_post_spike_ca_ip{:.2f}_N{:d}.pdf".format(ip3, N), transparent=True)
plt.show()