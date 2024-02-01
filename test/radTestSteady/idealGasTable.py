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


def T(R, e, rho):
    cv = 5/2 * R
    return e/cv


eCoeff = np.linspace(11.8748, 11.8748 + 0.1*40 , num=40)
print(eCoeff)
TCoeff = np.log(T(R, np.exp(eCoeff) , 0))

print(TCoeff)