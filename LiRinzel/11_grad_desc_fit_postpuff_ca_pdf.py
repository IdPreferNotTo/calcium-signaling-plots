import numpy as np
import matplotlib.pyplot as plt
from LiRinzel import utilities as ut


class LinIF:
    def __init__(self, ip3, r0, tau, xR, D):
        self.tau = tau
        self.xR = xR
        self.D = D

        self.r0 = r0
        self.x0 = ut.fixed_point(ip3)[0]
        self.x_max = ut.critical_point(ip3)[0]
        self.integral_caR_caT = self.integral_ca_caT(xR)
        self.var_integral_caR_caT = self.var_integral_ca_caT(xR)

    def update_par(self, r0, D, tau):
        self.D = D
        self.tau = tau
        self.r0 = r0
        self.integral_caR_caT = self.integral_ca_caT(xR)

    def U(self, x):
        return np.power(x - self.x0, 2) / (2 * self.tau)

    def integral_ca_caT(self, x):
        integral = 0
        dx = (self.x_max - x) / 100
        for x in np.linspace(x, self.x_max, 100):
            integral += dx * np.exp(self.U(x) / self.D)
        return integral

    def var_integral_ca_caT(self, x):
        integral = 0
        dx = (self.x_max - x) / 100
        for x in np.linspace(x, self.x_max, 100):
            integral += dx * (self.U(x) / self.D) * np.exp(self.U(x) / self.D)
        return integral

    def p0x(self, ca, r0, D, tau):
        integral = 0
        if ca < self.xR:
            dx = (self.x_max - self.xR) / 100
            for x in np.linspace(self.xR, self.x_max, 100):
                integral += dx * np.exp(np.power(x - self.x0, 2) / (2 * tau * D))
            p0 = (r0 / D) * np.exp(- np.power(ca - self.x0, 2) / (2 * tau * D)) * integral
        else:
            dx = (self.x_max - ca) / 100
            for x in np.linspace(ca, self.x_max, 100):
                integral += dx * np.exp(np.power(x - self.x0, 2) / (2 * tau * D))
            p0 = (r0 / D) * np.exp(-np.power(ca - self.x0, 2) / (2 * tau * D)) * integral
        return p0

    def p0(self, x):
        if x < self.xR:
            p0 = self.r0 / self.D * np.exp(-self.U(x) / self.D) * self.integral_caR_caT
        else:
            p0 = self.r0 / self.D * np.exp(-self.U(x) / self.D) * self.integral_ca_caT(x)
        return p0

    def gradient_descent(self, xs, ysaim):
        epsilon = 0.01
        sqerr = 0
        sqerr_dD = 0
        sqerr_dtau = 0
        sqerr_dr0 = 0
        m = len(xs)
        for x, yaim in zip(xs, ysaim):
            sqerr += (1 / 2) * (self.p0x(x, self.r0, self.D, self.tau) - yaim) ** 2
            sqerr_dD += (1 / 2) * (self.p0x(x, self.r0 , 1.01 * self.D, self.tau) - yaim) ** 2
            sqerr_dtau += (1 / 2) * (self.p0x(x, self.r0, self.D, 1.01 * self.tau) - yaim) ** 2
            sqerr_dr0 += (1 / 2) * (self.p0x(x, 1.01*self.r0, self.D, self.tau) - yaim) ** 2

        gradient_D = (1 / m)*(sqerr_dD - sqerr) / (0.01 * self.D)
        gradient_tau = (1 / m)*(sqerr_dtau - sqerr) / (0.01 * self.tau)
        gradient_r0 = (1 / m) * (sqerr_dr0 - sqerr) / (0.01 * self.r0)
        self.D -= epsilon * gradient_D
        self.tau -= epsilon * gradient_tau
        self.r0 -= epsilon * gradient_r0
        self.integral_caR_caT = self.integral_ca_caT(xR)


class ExpIF:
    def __init__(self, ip3, r0, tau, xR, xT, delta, D):
        self.tau = tau
        self.xR = xR
        self.xT = xT
        self.delta = delta
        self.D = D

        self.r0 = r0
        self.x0 = ut.fixed_point(ip3)[0]
        self.x_max = ut.critical_point(ip3)[0]

        self.integral_caR_caT = self.integral_ca_caT(xR)
        self.full_integral = self.integral_ca_caT(0)
        self.norm = self.normalization()


    def normalization(self):
        norm = 0
        xs = np.linspace(0, self.x_max, 200)
        dx = self.x_max/200
        for x in xs:
            p0 = self.p0(x)
            norm += p0*dx
        return norm


    def U(self, x):
        return np.power(x - self.x0, 2) / (2 * self.tau) - (self.delta ** 2 / self.tau) * np.exp(
            (x - self.xT) / self.delta)


    def integral_ca_caT(self, x):
        integral = 0
        dx = (self.x_max - x) / 100
        for x in np.linspace(x, self.x_max, 100):
            integral += dx * np.exp(self.U(x) / self.D)
        return integral


    def p0(self, x):
        if x < self.xR:
            p0 = self.r0 / self.D * np.exp(-self.U(x) / self.D) * self.integral_caR_caT
        else:
            p0 = self.r0 / self.D * np.exp(-self.U(x) / self.D) * self.integral_ca_caT(x)
        return p0


    def p0x(self, ca, xR, D, tau, xT, delta):
        integral = 0
        if ca < xR:
            integral_range = np.linspace(xR, self.x_max, 100, endpoint=True)
            dx = (self.x_max - xR) / 100
        else :
            integral_range = np.linspace(ca, self.x_max, 100, endpoint=True)
            dx = (self.x_max - ca) / 100

        for x in integral_range:
            integral += dx * np.exp(np.power(x - self.x0, 2) / (2 * tau * D) - (delta ** 2 / (tau * D)) * np.exp(
            (x - xT) / delta))
        p0 = (r0 / D) * np.exp(-np.power(ca - self.x0, 2) / (2 * tau * D) + (delta ** 2 / (tau * D)) * np.exp(
            (ca - xT) / delta)) * integral
        return p0


    def gradient_descent(self, xs, ysaim):
        epsilon = 0.01

        # Used to calculate the square error if a single variable is varied
        sqerr = 0
        sqerr_dxR = 0
        sqerr_dD = 0
        sqerr_dtau = 0
        sqerr_dxT = 0
        sqerr_ddelta = 0
        m = len(xs)

        # Ensures that the probability is normalized
        sum_p0xdxR = 1
        sum_p0xdD = 1
        sum_p0xdtau = 1
        sum_p0xdxT = 1
        sum_p0xddelta = 1

        dx = xs[1] - xs[0]
        for x in xs:
            sum_p0xdxR += self.p0x(x, 1.01*self.xR, self.D, self.tau, self.xT, self.delta)*dx
            sum_p0xdD += self.p0x(x, self.xR, 1.01*self.D, self.tau, self.xT, self.delta)*dx
            sum_p0xdtau += self.p0x(x, self.xR, self.D, 1.01*self.tau, self.xT, self.delta)*dx
            sum_p0xdxT += self.p0x(x, self.xR, self.D, self.tau, 1.01*self.xT, self.delta)*dx
            sum_p0xddelta += self.p0x(x, self.xR, self.D, self.tau, self.xT, 1.01*self.delta)*dx
        for x, yaim in zip(xs, ysaim):
            sqerr += (1 / 2) * (self.p0x(x, self.xR, self.D, self.tau, self.xT, self.delta) - yaim) ** 2
            sqerr_dxR += (1 / 2) * (self.p0x(x, 1.01*self.xR, self.D, self.tau, self.xT, self.delta)/sum_p0xdxR - yaim) ** 2
            sqerr_dD += (1 / 2) * (self.p0x(x, self.xR, 1.01*self.D, self.tau, self.xT, self.delta)/sum_p0xdD - yaim) ** 2
            sqerr_dtau += (1 / 2) * (self.p0x(x, self.xR, self.D, 1.01*self.tau, self.xT, self.delta)/sum_p0xdtau - yaim) ** 2
            sqerr_dxT += (1 / 2) * (self.p0x(x, self.xR, self.D, self.tau, 1.01*self.xT, self.delta)/sum_p0xdxT - yaim) ** 2
            sqerr_ddelta += (1 / 2) * (self.p0x(x, self.xR, self.D, self.tau, self.xT, 1.01*self.delta)/sum_p0xddelta - yaim) ** 2

        gradient_xR = (1 / m) * (sqerr_dxR - sqerr) / (0.01 * self.xR)
        gradient_D = (1 / m)*(sqerr_dD - sqerr) / (0.01 * self.D)
        gradient_tau = (1 / m)*(sqerr_dtau - sqerr) / (0.01 * self.tau)
        gradient_xT = (1 / m) * (sqerr_dxT - sqerr) / (0.01 * self.xT)
        gradient_delta = (1 / m) * (sqerr_ddelta - sqerr) / (0.01 * self.delta)

        self.xR -=epsilon * gradient_xR
        self.D -= epsilon * gradient_D
        self.tau -= epsilon * gradient_tau
        self.xT -= epsilon * gradient_xT
        self.delta -= epsilon * gradient_delta
        self.integral_caR_caT = self.integral_ca_caT(self.xR)
        norm = 0
        for x in xs:
            norm += dx*self.p0(x)
        self.norm = norm


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
            i += 1
        if i == len(starts):
            break
    return [t_puff_triggered_averages, ca_puff_triggered_averages]


def get_max(lists):
    maximum = 0
    for list in lists:
        if len(list) > 0:
            if maximum < max(list):
                maximum = max(list)
    return maximum


def get_min(lists):
    minimum = 1
    for list in lists:
        if len(list) > 0:
            if minimum > min(list):
                minimum = min(list)
    return minimum


def probability(valuess, min, max, bins=100):
    xs = [min + (i+1) * (max - min) / bins for i in range(bins)]
    pxs = [0 for _ in range(bins)]
    for values in valuess:
        for value in values:
            index = int((value - min) * bins / (max - min)) - 1
            pxs[index] += 1
    sum = np.sum(pxs)
    pxs = [x * (len(xs)/ (sum * (max - min))) for x in pxs]
    return xs, pxs


if __name__ == "__main__":
    ip3 = 0.3
    f, ax = plt.subplots(1, 1, figsize=(6, 9 / 2))
    bins = 100
    N = 64
    data = np.loadtxt(
        "/home/lukas/CLionProjects/li-rinzel-calcium-phd/out/lr_model_ca_waves_N{:d}_ip{:.2f}".format(N, ip3),
        skiprows=10_000)
    ts, cas, j1, n_open = np.transpose(data)
    ipis, starts, stops = ut.spike_times(ts, cas, n_open, ip3)
    ts_pre, cas_pre = pre_puff_ca(ts, cas, starts, stops)

    maximum = get_max(cas_pre)
    cas, pcas = probability(cas_pre, 0, maximum)

    #tau = 0.1
    #D = 0.1
    #Delta = 0.015
    #xR = 0.085
    #xT = 0.1
    #r0 = 1 / np.mean(ipis)

    #linIF = LinIF(ip3, r0, tau, xR, D)
    #p0 = []
    #for i in range(501):
    #    linIF.gradient_descent(cas, pcas)
    #    if i % 100 == 0:
    #        print("r0", linIF.r0, "tau: ", linIF.tau, "D: ", linIF.D)
    #        p0 = []
    #        sqerr = 0
    #        for x, pca in zip(cas, pcas):
    #            p0x = linIF.p0(x)
    #            p0.append(p0x)
    #        ax.plot(cas, p0)

    #r0 = linIF.r0
    #tau = linIF.tau
    #D = linIF.D

    #ax.plot(cas, pcas, label="r0={:.3f}, tau = {:.3f}, D= {:.3f}".format(r0, tau, D))
    #plt.legend()
    #plt.savefig("/home/lukas/Plots/lr_ca_phd/lr_linIF_grad_descent_ip{:.2f}_N{:d}.pdf".format(ip3, N), transparent=True)
    #plt.show()

    tau = 0.1
    D = 0.005
    Delta = 0.01
    xR = 0.09
    xT = 0.1
    r0 = 1 / np.mean(ipis)


    exIF = ExpIF(ip3, r0, tau, xR, xT, Delta, D)
    p0 = []
    norm = exIF.normalization()
    for x in cas:
       p0.append(exIF.p0(x)/norm)
    ax.plot(cas, p0)

    for i in range(0):
        exIF.gradient_descent(cas, pcas)
        print(i)
        if i % 20 == 0:
            print("tau: ", exIF.tau, "D: ", exIF.D, "xR:", exIF.xR, "xT: ", exIF.xT, "Delta: ", exIF.delta)
            p0 = []
            sqerr = 0
            for x, pca in zip(cas, pcas):
                p0x = exIF.p0(x)
                p0.append(p0x)
            ax.plot(cas, p0)

    ax.plot(cas, pcas, label=r"$r_0={:.3f}, \tau = {:.3f}, D={:.3f}, \Delta={:.3f}, x_R: {:.3f}, x_T={:.3f}$".format(exIF.r0, exIF.tau, exIF.D, exIF.delta, exIF.xR, exIF.xT), c="C3")
    plt.legend()
    plt.savefig("/home/lukas/Plots/lr_ca_phd/lr_expIF_grad_descent_ip{:.2f}_N{:d}.pdf".format(ip3, N), transparent=True)
    plt.show()