import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


l = 1                    # cm
c = 2.998e+8             # speed of light (m/s)
Tleft = 6029.62          # (K)
sigma = 5.670374419e-8   # W/m2/K4
a = 7.5657e-16           # J/cm3/K4
F = sigma * Tleft**4     # W/m2
absorptivity = 1         # m-1
alpha = 1e-15

print("time ",  alpha /sigma / absorptivity / 16 * 1 )
print("time ",  alpha /sigma / absorptivity / 16 * 100 )
print("tau : ", 4*a*c*1e-4*absorptivity/alpha*1e-6)
print("tau : ", 4*a*c*1e-4*absorptivity/alpha*1e-6)

tau = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100])
x = np.array([0, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 10, 15, 20])

#  u is the dimensionless radiation field, defined as: c/4 * G/F_inc, where F_inc is the incoming flux. Solution for tau = 100
u_100 = np.array([0.90895, 0.90107 , 0.88926, 0.86965, 0.85011, 0.83067, 0.71657, 0.54080, 0.38982, 0.26789, 0.10906, 0.03624 ])
v_100 = np.array([0.90849, 0.90057 , 0.88871, 0.86900, 0.84937, 0.82983, 0.71521, 0.53877, 0.38745, 0.26551, 0.10732, 0.03534 ])


#  u is the dimensionless radiation field, defined as: c/4 * G/F_inc, where F_inc is the incoming flux. Solution for tau = 1

u_0001 = np.array([0.03016, 0.00034])
u_0003 = np.array([0.05130, 0.00605, 0.00003])
u_001 = np.array([0.09040, 0.03241, 0.00361, 0.00001])
u_003 = np.array([0.14769, 0.08522, 0.03043, 0.00294, 0.00012])
u_01 = np.array([0.24023, 0.18003, 0.11024, 0.04111, 0.01217, 0.00280])
u_03 = np.array([0.34619, 0.29261, 0.22334, 0.13531, 0.07653, 0.04016, 0.00014])
u_1 = np.array([0.465099, 0.42133 , 0.36020, 0.27323, 0.20332, 0.14837, 0.01441, 0.00005, 0.00001])
u_100 = np.array([0.90895, 0.90107 , 0.88926, 0.86965, 0.85011, 0.83067, 0.71657, 0.54080, 0.38982, 0.26789, 0.10906, 0.03624 ])


v_0003 = np.array([0.00010, 0.00001])
v_001 = np.array([0.00062, 0.00014 , 0.000001])
v_003 = np.array([0.00302, 0.00135 , 0.00034, 0.00002])
v_01 = np.array([0.01641, 0.01068 , 0.00532, 0.00143, 0.00032, 0.00005])
v_03 = np.array([0.06844, 0.05353 , 0.03639, 0.01822, 0.00854, 0.00367, 0.00001])
v_1 = np.array([0.24762, 0.21614 , 0.17530 , 0.12182, 0.08306, 0.05556, 0.00324, 0.00001])
v_100 = np.array([0.90849, 0.90057 , 0.88871, 0.86900, 0.84937, 0.82983, 0.71521, 0.53877, 0.38745, 0.26551, 0.10732, 0.03534 ])



#  u is the dimensionless radiation field, defined as: c/4 * G/F_inc, where F_inc is the incoming flux. Solution for tau = 100



radiationCase_1em8 = np.loadtxt('postProcessing/sampleDict/1e-06/Centerline_G_T.xy')
radiationCase_1em6 = np.loadtxt('postProcessing/sampleDict/2e-06/Centerline_G_T.xy')
radiationCase_1em4 = np.loadtxt('postProcessing/sampleDict/3e-06/Centerline_G_T.xy')
# radiationCase_1em2 = np.loadtxt('postProcessing/sampleDict/1e-05/Centerline_G_T.xy')

z_1em8 = radiationCase_1em8[: , 0]
z_1em6 = radiationCase_1em6[: , 0]
z_1em4 = radiationCase_1em4[: , 0]
# z_1em2 = radiationCase_1em2[: , 0]

xNum_1em8 = np.sqrt(3) * absorptivity * z_1em8
xNum_1em6 = np.sqrt(3) * absorptivity * z_1em6
xNum_1em4 = np.sqrt(3) * absorptivity * z_1em4
# xNum_1em2 = np.sqrt(3) * absorptivity * z_1em2


G_1em8 = radiationCase_1em8[:, 1]      
G_1em6 = radiationCase_1em6[:, 1]      
G_1em4 = radiationCase_1em4[:, 1]      
# G_1em2 = radiationCase_1em2[:, 1]      

T_1em8 = radiationCase_1em8[:, 2]   
T_1em6 = radiationCase_1em6[:, 2]   



plt.plot(x[:2], u_0001 ,  label='u Analytical @ tau = 0.001')
plt.plot(x[:3], u_0003 ,  label='u Analytical @ tau = 0.003')
plt.plot(x[:4], u_001  ,  label='u Analytical @ tau = 0.01')
plt.plot(x[:5], u_003  ,  label='u Analytical @ tau = 0.03')
plt.plot(x[:6], u_01  ,  label='u Analytical @ tau = 0.1')
plt.plot(x[:7], u_03  ,  label='u Analytical @ tau = 0.3')
# plt.plot(x[:-3], u_1 ,  label='u Analytical @ tau = 1')
plt.plot(xNum_1em8   ,  G_1em8/4/F , 'x' , label='u Numerical @ t = 1e-6')
plt.plot(xNum_1em6  ,   G_1em6/4/F , 'x' , label='u Numerical @ t = 3e-6')
plt.plot(xNum_1em4  ,   G_1em4/4/F , 'x' , label='u Numerical @ t = 1e-5 ')
# plt.plot(xNum_1em2  ,   G_1em2/4/F , 'x' , label='u Numerical @ t = 1e-5')
plt.xlabel('Position (dimensionless)')
plt.ylabel('Radiation intensity (dimensionless)')
plt.grid()
plt.legend()
plt.show()


plt.plot(x[:-4], v_1, 'x' , label='u Analytical @ tau = 1')
plt.plot(xNum_1em8  ,  2 * sigma / c *(T_1em8)**4  , label='u Numerical @ tau = 1')
plt.xlabel('Position (dimensionless)')
plt.ylabel('Temperature (dimensionless)')

# plt.xscale('log')

plt.legend()
plt.grid()
plt.show()