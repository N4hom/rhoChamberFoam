import pandas as pd
import matplotlib.pyplot as plt

def timeIntegrate(time, data):
    
    integral = 0
    deltaT = time[0]
    for i in range(len(time)-1):
        integral = integral + data[i+1] * deltaT
        deltaT = time[i+1] - time[i]

    return integral

def mass(time, data):
    
    mass = []
    integral = 0
    deltaT = time[0]
    for i in range(len(time)-1):
        integral = integral + data[i+1] * deltaT
        deltaT = time[i+1] - time[i]
        mass.append(integral)

    return mass



# Read data into a pandas DataFrame
innerMassFluxData = pd.read_csv("postProcessing/innerMassFlux/0/surfaceFieldValue.dat", skiprows=6, sep='\s+')

timeInner = innerMassFluxData[innerMassFluxData.columns[0]]
innerMassFlux = innerMassFluxData[innerMassFluxData.columns[1]]

outerMassFluxData = pd.read_csv("postProcessing/outerMassFlux/0/surfaceFieldValue.dat", skiprows=6, sep='\s+')

timeOuter = outerMassFluxData[outerMassFluxData.columns[0]]
outerMassFlux = outerMassFluxData[outerMassFluxData.columns[1]]


print(timeIntegrate(timeInner , innerMassFlux))
print(timeIntegrate(timeOuter , outerMassFlux))

massTot = mass(timeOuter, outerMassFlux)


# Plot the data
#plt.plot(timeInner, innerMassFlux,  linestyle='-', color='b' , label='inner')
#plt.plot(timeOuter, outerMassFlux,  linestyle='-', color='r', label='outer')
plt.plot(timeOuter[1:], massTot,  linestyle='-', color='r', label='outer')
plt.legend()
plt.xlabel('time (s)')
plt.ylabel('Mass flux (kg/s)')
plt.grid(True)
plt.show()


