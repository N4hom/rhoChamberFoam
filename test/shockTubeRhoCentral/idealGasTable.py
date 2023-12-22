import numpy as np

def rho(R, p, T):
    return p / R / T

def e(R, p, T):
	return 3/2*R*1000*T

def cSqr(R , p ,T):
    return np.sqrt(5/3 * R * 1000 * T)

# Set the desired precision
precision = 3

pressureExp = np.linspace(10, 19, num=30)
print(pressureExp[1:] - pressureExp[:-1])
temperatureExp = np.linspace(5, 12, num=30)
print(temperatureExp[1:] - temperatureExp[:-1])

pressure = np.exp(pressureExp)
temperature = np.exp(temperatureExp)


densityTable = np.zeros((len(pressureExp), len(temperatureExp)))
energyTable = np.zeros((len(pressureExp), len(temperatureExp)))
cSqrTable = np.zeros((len(pressureExp), len(temperatureExp)))

for i in range(len(temperatureExp)):
    densityTable[i] = np.log(rho(28.9, pressure[i], temperature))
    energyTable[i] = np.log(e(28.9, pressure[i], temperature))
    cSqrTable[i] = np.log(cSqr(28.9 , pressure[i] , temperature))
    print(','.join(map(lambda x: f"{x:.{precision}f}", densityTable[i])))

print()

for i in range(len(temperatureExp)):
    print(', '.join(map(lambda x: f"{x:.{precision}f}", energyTable[i])))

print()

for i in range(len(temperatureExp)):
    print(', '.join(map(lambda x: f"{x:.{precision}f}", cSqrTable[i])))

print()
