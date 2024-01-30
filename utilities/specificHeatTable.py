import numpy as np

R_ = 8314.0
M = 28.9

R = R_/M

alpha = 1e-15

def rho(R, p, T):
    return p / R / T

def e(R, p, T):
	return 5/2*R*T

def cSqr(R , p ,T):
    return np.sqrt(5/3 * R * T)


def T(alpha, rho, e):
    return (4/alpha * e)**(1/4)

precision = 5


eLog = np.arange(-50, 5 , 1)

Tlog = np.log(T(alpha , 0, np.exp(eLog)))


print(', '.join(map(lambda x: f"{x:.{precision}f}", Tlog)))
   
print(np.exp(Tlog))

import matplotlib.pyplot as plt

plt.plot(np.exp(eLog) , np.exp(Tlog))
plt.show()