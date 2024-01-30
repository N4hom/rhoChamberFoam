import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


l = 1 # cm
c = 2.998e+8 # speed of light (m/s)
Tleft = 6029.62 # (K)
sigma = 5.670374419e-8   # W/m2/K4
F = sigma * Tleft**4    # W/m2
print("F " , F)
print((c/4/sigma)**(1/4))
absorptivity = 0.1 # m-1
alpha = 1e-15

time = alpha /sigma / absorptivity / 16 * 0.001
print(time)
alpha = 16 * sigma * absorptivity * time / 0.001   
print(alpha)

time = 1e-8
tau = 16 * sigma * absorptivity * time / alpha
print("tau " ,16 * sigma * absorptivity * time / alpha)
print("time ",  alpha /sigma / absorptivity / 16 * 1 )
print("time ",  alpha /sigma / absorptivity / 16 * 100 )

tau = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100])
x = np.array([0, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 10, 15, 20])

#  u is the dimensionless radiation field, defined as: c/4 * G/F_inc, where F_inc is the incoming flux. Solution for tau = 100
u_100 = np.array([0.90895, 0.90107 , 0.88926, 0.86965, 0.85011, 0.83067, 0.71657, 0.54080, 0.38982, 0.26789, 0.10906, 0.03624 ])
v_100 = np.array([0.90849, 0.90057 , 0.88871, 0.86900, 0.84937, 0.82983, 0.71521, 0.53877, 0.38745, 0.26551, 0.10732, 0.03534 ])

#  u is the dimensionless radiation field, defined as: c/4 * G/F_inc, where F_inc is the incoming flux. Solution for tau = 1
u_1 = np.array([0.465099, 0.42133 , 0.36020, 0.27323, 0.20332, 0.14837, 0.01441, 0.00005, 0.00001])
v_1 = np.array([0.24762, 0.21614 , 0.17530 , 0.12182, 0.08306, 0.05556, 0.00324, 0.00001])


# radiationCase_1em10 = np.loadtxt('postProcessing/sampleDict/1e-11/Centerline_G_T.xy')
# radiationCase_1em8 = np.loadtxt('postProcessing/sampleDict/1e-09/Centerline_G_T.xy')
radiationCase_1em8 = np.loadtxt('postProcessing/sampleDict/1.1e-08/Centerline_G_T.xy')
radiationCase_1em6 = np.loadtxt('postProcessing/sampleDict/1.111e-06/Centerline_G_T.xy')

# z_1em10 = radiationCase_1em10[: , 0]
z_1em8 = radiationCase_1em8[: , 0]
z_1em6 = radiationCase_1em6[: , 0]

# xNum_1em10 = np.sqrt(3) * absorptivity * z_1em10
xNum_1em8 = np.sqrt(3) * absorptivity * z_1em8
xNum_1em6 = np.sqrt(3) * absorptivity * z_1em6

# print(xNum_1em10)

# G_1em10 = radiationCase_1em10[:, 1]     # J/m3 --> J/cm3
G_1em8 = radiationCase_1em8[:, 1]       # J/m3 --> J/cm3
G_1em6 = radiationCase_1em6[:, 1]       # J/m3 --> J/cm3

# print("G_1em8 " , G_1em8[0])
# print(c/4 * G_1em8[0]/F)


# T_1em10 = radiationCase_1em10[:, 2] 
T_1em8 = radiationCase_1em8[:, 2]   
T_1em6 = radiationCase_1em6[:, 2]   



# uNum_1em10 = c/4 * G_1em10/F  
# uNum_1em8 = c/4 * G_1em8/F  
# uNum_1em8 = c/4 * G_1em8/F  



# plt.plot(xNum_1em10 , uNum_1em10 , label='t = 0.1 ns')
# plt.xscale('log')
# plt.plot(xNum_1em8 ,  uNum_1em8 , label='t = 1 ns')
# plt.xscale('log')
# plt.plot(xNum_1em8 ,  uNum_1em8 , label='t = 10 ns')
# plt.xscale('log')

# G_analytical = 4/c * F * u

plt.plot(x[:-3]/np.sqrt(3) / absorptivity, u_1, 'x' , label='u Analytical @ tau = 1')
plt.plot(x/np.sqrt(3) / absorptivity, u_100, 'x' , label='u Analytical @ tau = 100')
# plt.plot(x[:-3]/np.sqrt(3) / absorptivity, u_100, 'x' , label='u Analytical tau = 100')
plt.plot(z_1em8 * 10 ,  G_1em8/c , label='u Numerical @ tau = 1')
plt.plot(z_1em6 * 10 ,  G_1em6/c , label='u Numerical @ tau = 100')
plt.xlabel('Position(m)')
plt.ylabel('Radiation intensity (J/m3)')

# plt.xscale('log')

plt.legend()
plt.grid()
plt.show()

plt.plot(x[:-4]/np.sqrt(3) / absorptivity, v_1, 'x' , label='u Analytical @ tau = 1')
# plt.plot(x[:-3]/np.sqrt(3) / absorptivity, u_100, 'x' , label='u Analytical tau = 100')
plt.plot(z_1em8 * 10 ,  4 * sigma / c *(T_1em8)**4  , label='u Numerical @ tau = 1')
# plt.plot(z_1em6 ,  G_1em6/c , label='tau = 100')
plt.xlabel('Position(m)')
plt.ylabel('Radiation intensity (J/m3)')

# plt.xscale('log')

plt.legend()
plt.grid()
plt.show()