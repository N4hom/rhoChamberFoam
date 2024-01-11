import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Load the data from the file into an array
rhoCentralFoam = np.loadtxt('shockTubeRhoCentral/postProcessing/sampleDict/0.005/Centerline_T_p_rho_U.xy')
chamberFoam = np.loadtxt('myThermoTabTest/postProcessing/sampleDict/0.005/Centerline_T_p_rho_U.xy')


# Extract the first two columns into x and y arrays
axialDistanceRhoCentralFoam = rhoCentralFoam[:, 0]
axialDistanceChamberFoam = chamberFoam[:, 0]

temperatureRhoCentralFoam = rhoCentralFoam[: , 1]
temperatureChamberFoam = chamberFoam[: , 1]

pressureRhoCentralFoam = rhoCentralFoam[:, 2]
pressureChamberFoam = chamberFoam[:, 2]

densityRhoCentralFoam = rhoCentralFoam[:, 3]
densityChamberFoam = chamberFoam[:, 3]

velRhoCentralFoam = rhoCentralFoam[:, 4]
velChamberFoam = chamberFoam[:, 4]


# PRESSURE
plt.figure()
plt.plot(axialDistanceRhoCentralFoam, pressureRhoCentralFoam , label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam, pressureChamberFoam , label='chamberFoam')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('Pressure (Pa)', fontsize=13)
plt.legend()
plt.grid()
plt.show()

# TEMPERATURE
plt.figure()
plt.plot(axialDistanceRhoCentralFoam, temperatureRhoCentralFoam, label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam, temperatureChamberFoam , label='chamberFoam')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('Temperature (K)', fontsize=13)
plt.legend()
plt.grid()
plt.show()

# DENSITY
plt.figure()
plt.plot(axialDistanceRhoCentralFoam, densityRhoCentralFoam, label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam, densityChamberFoam , label='chamberFoam')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('Density (kg/m3) ' , fontsize=13)
plt.legend()
plt.grid()
plt.show()

# VELOCITY
plt.figure()
plt.plot(axialDistanceRhoCentralFoam, velRhoCentralFoam, label='rhoCentralFoam')
plt.plot(axialDistanceChamberFoam, velChamberFoam , label='chamberFoam')
plt.xlabel('x (a.u.)', fontsize=13)
plt.ylabel('Velocity (m/s) ' , fontsize=13)
plt.legend()
plt.grid()
plt.show()
