import numpy as np

R_ = 8314.0   #J/K/kmol
gamma = 1.67
M = 7.2947

R = R_/M   #J/kg/K
precision = 9


def rho(R, p, T):
    return p / R / T

def e(R, p, T):
	return 5/2*R*T

def cSqr(R , p ,T):
    return np.sqrt(5/3 * R * T)

def Cv(rho, e, R):
    return 1.5*R

def p(rho, e):
    return (gamma - 1) * rho * e

def T(rho, e, R):
    return e/Cv(e,rho,R)


rhoLog = np.linspace(-5, 3, num=50)

eLog = np.linspace(12, 23, 50)
print(eLog[1:]-eLog[:-1])
print(rhoLog[1:]-rhoLog[:-1])



print("temperature coefficients: ")
for i in range(len(rhoLog)):
    print(';'.join(map(lambda x: f"{x:.{precision}f}", np.log(T(10**(rhoLog[i]), np.exp(eLog), R)))))
print()

print("Pressure coefficients: ")
for i in range(len(rhoLog)):
    print(';'.join(map(lambda x: f"{x:.{precision}f}", np.log(p(10**(rhoLog[i]), np.exp(eLog))))))
print()

print("Pressure values: ")
for i in range(len(rhoLog)):
    print(';'.join(map(lambda x: f"{x:.{precision}f}", p(10**(rhoLog[i]), np.exp(eLog)))))
print()

print("Temperature values: ")
for i in range(len(rhoLog)):
    print(';'.join(map(lambda x: f"{x:.{precision}f}", T(10**(rhoLog[i]), np.exp(eLog), R))))
print()