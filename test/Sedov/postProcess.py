import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

dr = 0.005
gamma = 1.67
p0 = 1e5
rho0 = 1.16
p1 = 6.54822e+10
E = p1/rho0/(gamma - 1) * 4/3 * np.pi * dr**3
A = np.sqrt(E/5.36/rho0)
a = 344.0

print(A)
# rhoCentralFoam
rhoCentralFoam = np.loadtxt('SedovTest/postProcessing/sampleDict/0.0001/Centerline_p_rho_U.xy')

axialDistanceRhoCentralFoam = rhoCentralFoam[:, 0]
pressureRhoCentralFoam = rhoCentralFoam[:, 1]
densityRhoCentralFoam = rhoCentralFoam[:, 2]
velRhoCentralFoam = rhoCentralFoam[:, 3]

# chamberFoam
rhoChamberFoamHLLC = np.loadtxt('SedovHLLC/postProcessing/sampleDict/0.0001/Centerline_p_rho_U.xy')

axialDistanceChamberFoam = rhoChamberFoamHLLC[:, 0]
pressureChamberFoam = rhoChamberFoamHLLC[:, 1]
densityChamberFoam = rhoChamberFoamHLLC[:, 2]
velChamberFoam = rhoChamberFoamHLLC[:, 3]

blastFoam = np.loadtxt('Sedov_3D/postProcessing/sampleDict/0.0001/Centerline_p_rho_U.xy')

axialDistanceBlastFoam = blastFoam[:, 0]
pressureBlastFoam = blastFoam[:, 1]
densityBlastFoam = blastFoam[:, 2]
velBlastFoam = blastFoam[:, 3]



f = np.array([1.167 , 0.949 , 0.808 , 0.711 , 0.643 ,  0.593  ,  0.556 ,  0.528  ,  0.507 ,  0.491 ,  0.478 ,  0.468 ,  0.461 , 0.455 ,  0.450 ,  0.447 , 0.444 , 0.442 , 0.440 , 0.439 , 0.438 , 0.438 , 0.437 ,  0.437 , 0.437 , 0.436 , 0.436])  # Adding one data point just to make the plot better
f = np.flip(f)

eta = np.array([1.00 , 0.98 , 0.96 , 0.94 , 0.92 ,  0.90  ,  0.88 ,  0.86  ,  0.84 ,  0.82,  0.80 ,  0.78 ,  0.76 , 0.74 ,  0.72 ,  0.70 , 0.68 , 0.66 , 0.64 , 0.62 , 0.60 , 0.58 , 0.56 ,  0.54 , 0.52 , 0.50 , 0])
eta = np.flip(eta)

psi = np.array([0.007 ,0.007, 0.010 , 0.014 , 0.019 , 0.026 , 0.034 , 0.044 , 0.058 , 0.074 , 0.095 , 0.120 , 0.152 , 0.191 , 0.239 , 0.297 , 0.370 , 0.462 , 0.578 , 0.727 , 0.919 , 1.177 , 1.534 , 2.052 , 2.808 , 4.0 , 6.0])

eta_166 = np.array([0.0 , 0.50 , 0.7 ,  0.80 , 0.90 , 0.95 , 1.00])
psi_166 = np.array([0.0 , 0.05 , 0.29 , 0.63 , 1.14 , 2.30 , 4.00])

idxMaxP = np.argmax(pressureChamberFoam)

print("Pressure at the origin (5 cells):")
print("blastFoam : ", pressureBlastFoam[0]/pressureBlastFoam[idxMaxP])
print("rhoCentralFoam : ", pressureRhoCentralFoam[0]/pressureRhoCentralFoam[idxMaxP])
print("chamberFoam : ", pressureChamberFoam[0]/pressureChamberFoam[idxMaxP])
print("Analytical : " , f[0]/max(f))

plt.figure(1)
plt.plot(axialDistanceBlastFoam/axialDistanceBlastFoam[idxMaxP], (pressureBlastFoam)/pressureBlastFoam[idxMaxP] , label='blastFoam')
plt.plot(axialDistanceRhoCentralFoam/axialDistanceRhoCentralFoam[idxMaxP], (pressureRhoCentralFoam)/max(pressureRhoCentralFoam)  , label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam/axialDistanceChamberFoam[idxMaxP], (pressureChamberFoam)/max(pressureChamberFoam) , label='rhoChamberFoam HLLC')
plt.plot(eta, f/max(f) , label='analytical')
plt.xlim([0,1])
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('P/P0 ', fontsize=13)
plt.legend()
plt.grid()
plt.show()





plt.figure(1)
plt.plot(axialDistanceBlastFoam/0.25, pressureBlastFoam/max(pressureBlastFoam) * 1.167 , label='blastFoam')
plt.plot(axialDistanceRhoCentralFoam/0.25, pressureRhoCentralFoam/max(pressureRhoCentralFoam) * 1.167 , label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam/0.25, pressureChamberFoam/max(pressureChamberFoam) * 1.167, label='chamberFoam')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('P/Pmax ', fontsize=13)
plt.legend()
plt.grid()
plt.show()

# p/pmax * R^3 is supposed to be f1 

plt.figure(1)
plt.plot(axialDistanceBlastFoam/0.25, pressureBlastFoam/pressureBlastFoam[-1] , label='blastFoam')
plt.plot(axialDistanceRhoCentralFoam/0.25, pressureRhoCentralFoam/pressureRhoCentralFoam[-1] , label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam/0.25, pressureChamberFoam/pressureChamberFoam[-1] , label='chamberFoam')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('P/Pmax ', fontsize=13)
plt.legend()
plt.grid()
plt.show()



plt.figure(1)
plt.plot(axialDistanceBlastFoam/0.25, densityBlastFoam/rho0 , label='blastFoam')
plt.plot(axialDistanceRhoCentralFoam/0.25, densityRhoCentralFoam/rho0  , label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam/0.25, densityChamberFoam/rho0 , label='chamberFoam')
plt.plot(eta_166, psi_166 , 'o' , label='Taylor')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('rho (kg/m3)', fontsize=13)
plt.legend()
plt.grid()
plt.show()


plt.figure(1)
plt.plot(axialDistanceBlastFoam, velBlastFoam , label='blastFoam')
plt.plot(axialDistanceRhoCentralFoam, velRhoCentralFoam  , label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam, velChamberFoam , label='chamberFoam')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('velocity (m/s)', fontsize=13)
plt.legend()
plt.grid()
plt.show()




# rhoCentralFoam with only 1 cell as high pressure region
rhoCentralFoam1cell = np.loadtxt('SedovTestRhoCentral1cell/postProcessing/sampleDict/1e-05/Centerline_p_rho_U.xy')

axialDistanceRhoCentralFoam1cell = rhoCentralFoam1cell[:, 0]
pressureRhoCentralFoam1cell = rhoCentralFoam1cell[:, 1]
densityRhoCentralFoam1cell = rhoCentralFoam1cell[:, 2]
velRhoCentralFoam1cell = rhoCentralFoam1cell[:, 3]
idxMaxP1cell = np.argmax(pressureRhoCentralFoam1cell)

plt.figure(1)
plt.plot(axialDistanceRhoCentralFoam/axialDistanceRhoCentralFoam[idxMaxP], (pressureRhoCentralFoam)/max(pressureRhoCentralFoam)  , label='rhoCentralFoam 5 cells')
plt.plot(axialDistanceRhoCentralFoam1cell/axialDistanceRhoCentralFoam1cell[idxMaxP1cell], (pressureRhoCentralFoam1cell)/max(pressureRhoCentralFoam1cell)  , label='rhoCentralFoam 1 cell')
plt.xlim([0,1])
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('P/P0 ', fontsize=13)
plt.legend()
plt.grid()
plt.show()



plt.figure(1)
plt.plot(axialDistanceRhoCentralFoam1cell/axialDistanceRhoCentralFoam1cell[idxMaxP1cell], (pressureRhoCentralFoam1cell)/max(pressureRhoCentralFoam1cell) , 'x' , label='rhoCentralFoam')
plt.plot(eta, f/max(f) , label='Analytical')
plt.xlim([0,1])
plt.xlabel('x (a.u.)', fontsize=16)
plt.ylabel('$P/P_{max}$ ', fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.show()


plt.figure(1)
plt.plot(axialDistanceRhoCentralFoam1cell/axialDistanceRhoCentralFoam1cell[idxMaxP1cell], (densityRhoCentralFoam1cell)/max(densityRhoCentralFoam1cell) , 'x' , label='rhoCentralFoam')
plt.plot(eta, psi/max(psi)  ,label='Analytical')
plt.xlim([0,1])
plt.xlabel('x (a.u.)', fontsize=16)
plt.ylabel(r'$\rho/\rho_{max}$', fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.show()

plt.figure(1)
plt.plot(axialDistanceRhoCentralFoam1cell/axialDistanceRhoCentralFoam1cell[idxMaxP1cell], (densityRhoCentralFoam1cell)/densityRhoCentralFoam1cell[-1] , 'x' , label='rhoCentralFoam')
plt.plot(eta, psi  ,label='Analytical')
plt.xlim([0,1])
plt.xlabel('x (a.u.)', fontsize=16)
plt.ylabel(r'$\rho/\rho_0$', fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.show()



print("Pressure at the origin (1 cell):")
print("rhoCentralFoam : ", pressureRhoCentralFoam1cell[0]/pressureRhoCentralFoam1cell[idxMaxP1cell])
print("Analytical : " , f[0]/max(f))


# plt.figure(1)
# plt.xlabel('x (a.u.)', fontsize=13)
# plt.ylabel('P/P0', fontsize=13)
# plt.legend()
# plt.grid()
# plt.show()






# # PRESSURE
# plt.figure()
# plt.plot(axialDistanceRhoCentralFoam, pressureRhoCentralFoam , label='rhoCentralFoam')
# plt.plot(axialDistanceChamberFoam, pressureChamberFoam , label='chamberFoam')
# plt.xlabel('x (a.u.)', fontsize=13)
# plt.ylabel('Pressure (Pa)', fontsize=13)
# plt.legend()
# plt.grid()
# plt.show()

# # TEMPERATURE
# plt.figure()
# plt.plot(axialDistanceRhoCentralFoam, temperatureRhoCentralFoam, label='rhoCentralFoam')
# plt.plot(axialDistanceChamberFoam, temperatureChamberFoam , label='chamberFoam')
# plt.xlabel('x (a.u.)', fontsize=13)
# plt.ylabel('Temperature (K)', fontsize=13)
# plt.legend()
# plt.grid()
# plt.show()

# # DENSITY
# plt.figure()
# plt.plot(axialDistanceRhoCentralFoam, densityRhoCentralFoam, label='rhoCentralFoam')
# plt.plot(axialDistanceChamberFoam, densityChamberFoam , label='chamberFoam')
# plt.xlabel('x (a.u.)', fontsize=13)
# plt.ylabel('Density (kg/m3) ' , fontsize=13)
# plt.legend()
# plt.grid()
# plt.show()

# # VELOCITY
# plt.figure()
# plt.plot(axialDistanceRhoCentralFoam, velRhoCentralFoam, label='rhoCentralFoam')
# plt.plot(axialDistanceChamberFoam, velChamberFoam , label='chamberFoam')
# plt.xlabel('x (a.u.)', fontsize=13)
# plt.ylabel('Velocity (m/s) ' , fontsize=13)
# plt.legend()
# plt.grid()
# plt.show()
