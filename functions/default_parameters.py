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


def load_traces_markov(tau, p, cer, taua = 0, ampa=0, ip3 = 1.0, K=10, N=5):
    folder_markov = home + f"/Data/calcium_spikes_markov/"
    if not cer:
        folder_adap = "Data_no_adap/"
        file_traces = f"ca_markov_ip{ip3:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    else:
        folder_adap = f"Data_adap/tau{tau:.2e}_j{p:.2e}/"
        file_traces = f"ca_markov_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_markov + folder_adap + file_traces)


def load_traces_fixed_ci_markov(tau, p, ca, ip3 = 1.0, K=10, N=5):
    folder_markov_no_adap = home + f"/Data/calcium_spikes_markov/ca_fix/"
    file_traces = f"ca_markov_cafix{ca:.2f}_ip{ip3:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_markov_no_adap + file_traces)


def load_spike_times_markov(tau, p, cer, taua = 0, ampa=0, ip3 = 1.0, K=10, N=5):
    folder_markov = home + f"/Data/calcium_spikes_markov/"
    if not cer:
        folder_adap = "Data_no_adap/"
        file_spike_times = f"spike_times_markov_ip{ip3:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    else:
        folder_adap = f"Data_adap/tau{tau:.2e}_j{p:.2e}/"
        file_spike_times = f"spike_times_markov_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_markov + folder_adap + file_spike_times)


def load_spike_times_langevin(tau, p, cer, interpretation="strat", taua = 0, ampa=0, ip3 = 1.0, K=10, N=5):
    folder_markov = home + f"/Data/calcium_spikes_langevin_{interpretation}/"
    if not cer:
        folder_adap = "Data_no_adap/"
        file_spike_times = f"spike_times_langevin_ip{ip3:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}_0.dat"
    else:
        folder_adap = f"Data_adap/tau{tau:.2e}_j{p:.2e}/"
        file_spike_times = f"spike_times_langevin_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_markov + folder_adap + file_spike_times)

def load_ci_density(tau, p, s=1, K=10):
    folder = home + f"/Data/calcium_spikes_theory/probability_density/"
    file = f"p0_ip{s:.2f}_tau{tau:.2e}_j{p:.2e}_K{K:d}.dat"
    return np.loadtxt(folder + file)


def load_spike_times_markov_transient(tau, p, taua, ampa, ip3=1.0, K=10, N=5):
    folder_markov = home + f"/Data/calcium_spikes_markov/"
    folder_adap = f"Data_adap_transient/"
    file_spike_times = f"transient_spike_times_markov_ip{ip3:.2f}_taua{taua:.2e}_ampa{ampa:.2e}_tau{tau:.2e}_j{p:.2e}_K{K:d}_{N:d}.dat"
    return np.loadtxt(folder_markov + folder_adap + file_spike_times)