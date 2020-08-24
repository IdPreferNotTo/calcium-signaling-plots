import numpy as np
from LiRinzel import utilities as ut

v1 = 6
v2 = 0.11
v3 = 0.9
c0 = 2.
c1 = 0.185

k3 = 0.1
a2 = 0.2
d1 = 0.13
d2 = 1.049
d3 = 0.9434
d5 = 0.08234

ip3 = 0.3
n= ip3/(ip3+d1)
ca0, h0 = ut.fixed_point(ip3)
print(ca0, h0)
Q2 = d2*(ip3+d1)/(ip3 + d3)
J = np.array([[-v1*(n**3)*(h0**3)*(ca0**2)*(ca0*(4*d5 + ca0)*(c1 +1) -3*d5*c0)/((ca0+d5)**4) -v2*(c1+1) -2*v3*(k3**2)*ca0/((k3**2 +ca0**2)**2), -3*v1*(n**3)*((ca0/(ca0 +d5))**3)*(h0**2)*(ca0*(c1 +1) - c0)], [-a2*h0, -a2*(ca0+Q2)]])
w, v = np.linalg.eigvals(J)
print(w, v)