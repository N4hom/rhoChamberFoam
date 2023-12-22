import numpy as np

R_ = 8314.0
M = 28.9

R = R_/M

def rho(R, p, T):
    return p / R / T

def e(R, p, T):
	return 5/2*R*T

def cSqr(R , p ,T):
    return np.sqrt(5/3 * R * T)

scale = 'lin'

# Set the desired precision
precision = 3


if scale == 'lin':
    #pressure = np.linspace(np.exp(9.1), np.exp(11.52), num=10)
    pressure = np.array([1e5, 100100, 100200, 100300])
    print("Delta pressure ",pressure[1:] - pressure[:-1])
    temperature = np.linspace(np.exp(5), np.exp(6), num=10)
    print("Delta temperature ", temperature[1:] - temperature[:-1])
    

if scale == 'log':
    pressureExp = np.linspace(9, 13, num=4)
    print(pressureExp[1:] - pressureExp[:-1])
    temperatureExp = np.linspace(5, 12, num=4)
    print(temperatureExp[1:] - temperatureExp[:-1])
    pressure = np.exp(pressureExp)
    temperature = np.exp(temperatureExp)




densityTable = np.zeros((len(pressure), len(temperature)))
energyTable = np.zeros((len(pressure), len(temperature)))
cSqrTable = np.zeros((len(pressure), len(temperature)))


if scale == 'log':
    for i in range(len(temperature)):
        densityTable[i] = np.log(rho(R, pressure[i], temperature))
        energyTable[i] = np.log(e(R, pressure[i], temperature))
        cSqrTable[i] = np.log(cSqr(R , pressure[i] , temperature))

if scale == 'lin':
    for i in range(len(temperature)):
        densityTable[i] = rho(R, pressure[i], temperature)
        energyTable[i] = e(R, pressure[i], temperature)
        cSqrTable[i] = cSqr(R , pressure[i] , temperature)
    

print("Pressure scale : ", pressure)
print()
print("Temperature scale : ", temperature)

print("Density table coefficients:")
for i in range(len(temperature)):
    print(', '.join(map(lambda x: f"{x:.{precision}f}", densityTable[i])))
print()

if scale == 'log':
    print("Density table :")
    for i in range(len(temperature)):
        print(', '.join(map(lambda x: f"{x:.{precision}f}", np.exp(densityTable[i]))))
    print()
    print("Energy table :")
    for i in range(len(temperature)):
        print(', '.join(map(lambda x: f"{x:.{precision}f}", np.exp(energyTable[i]))))
    print()
    print("Speed of sound table: ")
    for i in range(len(temperature)):
        print(', '.join(map(lambda x: f"{x:.{precision}f}", np.exp(cSqrTable[i]))))

print()
print("Energy table coefficients:")
for i in range(len(temperature)):
    print(', '.join(map(lambda x: f"{x:.{precision}f}", energyTable[i])))

print()

print("Speed of sound table coefficients: ")
for i in range(len(temperature)):
    print(', '.join(map(lambda x: f"{x:.{precision}f}", cSqrTable[i])))
print()
