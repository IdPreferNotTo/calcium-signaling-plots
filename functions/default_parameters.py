import numpy as np
import os

home = os.path.expanduser("~")

r_ref = 20.
r_cls = 50.

def r_opn(ci, s=1, N=5):
    cR = 0.2
    cR3 = np.power(cR, 3)
    s3 = np.power(s, 3)
    r_opn_hat = 0.1 * (1  + cR3)/cR3 * (1 + s3)/s3
    r_opn = N * r_opn_hat * np.power(ci, 3)/(1 + np.power(ci, 3)) * np.power(s, 3)/(1 + np.power(s, 3))
    return r_opn

def load_traces_fixed_ci_markov(tau, p, ca, ip3 = 1.0, K=10, N=5):
    folder_markov_no_adap = home + f"/Data/calcium/markov/static/"
    file_traces = f"ca_markov_cafix{ca:.2f}_ip{ip3:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_markov_no_adap + file_traces)

def load_traces_markov(tau, p, cer=False, taua = 0, ampa=0, ip3 = 1.0, K=10, N=5):
    folder_markov = home + f"/Data/calcium/markov/"
    if not cer:
        folder = "renewal/traces/"
        file = f"ca_markov_ip{ip3:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    else:
        folder = f"adaptive/traces/"
        file = f"ca_markov_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_markov + folder + file)

def load_spike_times_markov(tau, p, cer, taua = 0, ampa=0, ip3 = 1.0, K=10, N=5):
    folder_markov = home + f"/Data/calcium/markov/"
    if not cer:
        folder = "renewal/sequence/"
        file = f"spike_times_markov_ip{ip3:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    else:
        folder = f"adaptive/sequence/"
        file = f"spike_times_markov_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    print(folder_markov + folder + file)
    return np.loadtxt(folder_markov + folder + file)


def load_spike_times_langevin(tau, p, cer, taua = 0, ampa=0, ip3 = 1.0, K=10, N=5, interpretation="strat",):
    folder_langevin= home + f"/Data/calcium/langevin/"
    if not cer:
        folder = "renewal/sequence/"
        file = f"spike_times_langevin_ip{ip3:.2f}_tau{tau:.2e}_j{p:.3e}_K{K:d}_5.dat"
    else:
        folder = f"adaptive/sequence/"
        file = f"spike_times_langevin_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_langevin + folder + file)

def load_adaptation_markov_transient(tau, p, taua, ampa, ip3=1.0, K=10, N=5):
    folder = home + f"/Data/calcium/markov/adaptive/transient/"
    file = f"transient_adaptation_markov_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder + file)

def load_spike_times_markov_transient(tau, p, taua, ampa, ip3=1.0, K=10, N=5):
    folder = home + f"/Data/calcium/markov/adaptive/transient/"
    file = f"transient_spike_times_markov_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder + file)

def load_spike_times_langevin_transient(tau, p, taua, ampa, ip3=1.0, K=10, N=5):
    folder = home + f"/Data/calcium/langevin/adaptive/transient/"
    file = f"transient_spike_times_langevin_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder + file)

def load_ci_density(tau, p, s=1, K=10):
    folder = home + f"/Data/calcium/theory/probability_density/"
    file = f"p0_ip{s:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}.dat"
    return np.loadtxt(folder + file)