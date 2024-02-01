import numpy as np 
import matplotlib.pyplot as plt
gamma = 1.4
rho0 = 1
p0 = 1e5
pExp = 6.54822e+10 # same value as in 0 folder
theta = np.linspace(0,1,num=100)
t = 1e-4
E = pExp/(gamma - 1)/rho0 # [J]
a = np.sqrt(E/5.36/rho0)

x = a * t**(2/5) * theta**((gamma - 1)/(2*gamma + 1)) * ((theta + 1)/2)**(-2/5) * ((theta + gamma)/(gamma + 1))**((gamma + 1)/(3*gamma - 1)) * ((3*(2 - gamma)*gamma + (2*gamma + 1))/(7 - gamma))**((-13*gamma**2 - 7*gamma + 12)/(5*(2*gamma + 1)*(3*gamma - 1)))

rho = (gamma + 1)/(gamma - 1) * rho0 * theta**(3/(2*gamma + 1)) * ((theta + gamma)/(gamma + 1)) * (-4/(3*gamma - 1)) * ((3 * (2 - gamma) * theta + (2*gamma + 1))/(7 - gamma))**((13 * gamma**2 - 7*gamma + 12)/((2 - gamma)*(2*gamma + 1)*(3 * gamma - 1)))

p = 8/(25 * (gamma + 1)) * rho0 * a**2 * t**-(6/5) * ((theta + 1)/2) ** 6/5 * ((theta + gamma)/(gamma + 1)) ** -4*gamma/(3*gamma - 1) * ((3 * (2 - gamma) * theta + (2*gamma + 1))/(7 - gamma)) ** ((13 * gamma **2 - 7 * gamma + 12)/(5*(2 - gamma) * (3 * gamma - 1)))



R = a * t **(2/5) * E ** 1/5 * rho0 **(-1/5)
f = np.array([1.167 , 0.949 , 0.808 , 0.711 , 0.643 ,  0.593  ,  0.556 ,  0.528  ,  0.507 ,  0.491 ,  0.478 ,  0.468 ,  0.461 , 0.455 ,  0.450 ,  0.447 , 0.444 , 0.442 , 0.440 , 0.439 , 0.438 , 0.438 , 0.437 ,  0.437 , 0.437 , 0.436 ])
f = np.flip(f)

eta = np.array([1.00 , 0.98 , 96 , 0.94 , 0.92 ,  0.90  ,  0.88 ,  0.86  ,  0.84 ,  0.82,  0.80 ,  0.78 ,  0.76 , 0.74 ,  0.72 ,  0.70 , 0.68 , 0.66 , 0.64 , 0.62 , 0.60 , 0.58 , 0.56 ,  0.54 , 0.52 , 0.50 ])
eta = np.flip(f)

plt.plot(f * R**-3 * p0 )
plt.grid()
plt.show()