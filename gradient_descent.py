import numpy as np
import matplotlib.pyplot as plt

class Target():
    def __init__(self, a, b, x0, x1):
        self.a = a
        self.b = b
        self.x0 = x0
        self.x1 = x1

    def function(self, x):
        return self.a*(x-self.x0) + self.b*(x-self.x1)**2

    def dfda(self, x):
        return (x-self.x0)

    def dfdb(self, x):
        return (x-self.x1)**2

    def dfdx0(self, x):
        return -self.a

    def dfdx1(self, x):
        return -2 *self.b * self.x1

    def update_par(self, da, db, dx0, dx1):
        self.a += da
        self.b += db
        self.x0 += dx0
        self.x1 += dx1


class GradientDesecent():
    def __init__(self, model):
        self.epsilon = 0.01
        self.model = model

    def opt(self, xs, target):
        da = 0
        db = 0
        dx0 = 0
        dx1 = 0
        cost = 0
        m = len(xs)
        for x, y in zip(xs, target):
            dif = (self.model.function(x) - y)/m
            da += dif * self.model.dfda(x)
            db += dif * self.model.dfdb(x)
            dx0 += dif * self.model.dfdx0(x)
            dx1 += dif * self.model.dfdx1(x)

        da = -self.epsilon*da/m
        db = -self.epsilon*db/m
        dx0 = -self.epsilon*dx0/m
        dx1 = -self.epsilon*dx1/m
        self.model.update_par(da, db, dx0, dx1)


if __name__ == "__main__":

    a = 2.
    b = 1
    x0 = 5.
    x1 = -1
    targetfunc = Target(a, b, x0, x1)
    xs = np.linspace(-10, 10, 200)
    ys = [targetfunc.function(x) for x in xs]

    plt.plot(xs, ys, c="C3")
    a_init = 5
    b_init = 2
    x0_init = 5
    x1_init = -1
    optfunc = Target(a_init, b_init, x0_init, x1_init)
    gradientdes = GradientDesecent(optfunc)
    optys = [optfunc.function(x) for x in xs]
    plt.plot(xs, optys)
    for i in range(10000):
        gradientdes.opt(xs, ys)
        if i%2000==0:
            optys = [optfunc.function(x) for x in xs]
            print(optfunc.a, optfunc.b, optfunc.x0, optfunc.x1)
            plt.plot(xs, optys)
    plt.show()



