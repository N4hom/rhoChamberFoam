import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# DATA FOR LIQUID FLIBE AND GEOMETRY
rho = 2000
H = 1
Area = 0.006282753 #m2

def calculateVelocity(forceData , timeData , rho , H , Area):
	impulse = []
	velocity = []
	impulseI = 0
	deltaT = timeData[0]
	for i in range(1,len(timeData)):
		deltaT = timeData[i] - timeData[i-1]
		impulseI = impulseI + forceData[i]*deltaT
		impulse.append(impulseI)
		velocity.append(impulseI/rho/H/Area)

	return velocity, impulse

# Read the data from the file
forceData = pd.read_csv('force.dat', delimiter='; ')

# Extract the columns
time = forceData['time(s)'].to_list()

forceSlabsx = np.array(forceData['slabs.x'].to_list())
forceSlabsy = np.array(forceData['slabs.y'].to_list())
forceSlabsz = np.array(forceData['slabs.z'].to_list())



forceJetsx1 = np.array(forceData['row_1.x'].to_list())
forceJetsy1 = np.array(forceData['row_1.y'].to_list())
forceJetsz1 = np.array(forceData['row_1.z'].to_list())


forceJetsx2 = np.array(forceData['row_2.x'].to_list())
forceJetsy2 = np.array(forceData['row_2.y'].to_list())
forceJetsz2 = np.array(forceData['row_2.z'].to_list())



forceJetsx3 = np.array(forceData['row_3.x'].to_list())
forceJetsy3 = np.array(forceData['row_3.y'].to_list())
forceJetsz3 = np.array(forceData['row_3.z'].to_list())


forceJetsx4 = np.array(forceData['row_4.x'].to_list())
forceJetsy4 = np.array(forceData['row_4.y'].to_list())
forceJetsz4 = np.array(forceData['row_4.z'].to_list())



forceJetsx5 = np.array(forceData['row_5.x'].to_list())
forceJetsy5 = np.array(forceData['row_5.y'].to_list())
forceJetsz5 = np.array(forceData['row_5.z'].to_list())


forceJetsx6 = np.array(forceData['row_6.x'].to_list())
forceJetsy6 = np.array(forceData['row_6.y'].to_list())
forceJetsz6 = np.array(forceData['row_6.z'].to_list())




forceJetsx7 = np.array(forceData['row_7.x'].to_list())
forceJetsy7 = np.array(forceData['row_7.y'].to_list())
forceJetsz7 = np.array(forceData['row_7.z'].to_list())



# Negligible but it could be useful for future 3D simulations


forceTotJets1 = np.sqrt(forceJetsx1**2 + forceJetsy1**2 + forceJetsz1**2)

forceTotJets2 = np.sqrt(forceJetsx2**2 + forceJetsy2**2 + forceJetsz2**2)

forceTotJets3 = np.sqrt(forceJetsx3**2 + forceJetsy3**2 + forceJetsz3**2)

forceTotJets4 = np.sqrt(forceJetsx4**2 + forceJetsy4**2 + forceJetsz4**2)

forceTotJets5 = np.sqrt(forceJetsx5**2 + forceJetsy5**2 + forceJetsz5**2)

forceTotJets6 = np.sqrt(forceJetsx6**2 + forceJetsy6**2 + forceJetsz6**2)


forceTotJets7 = np.sqrt(forceJetsx7**2 + forceJetsy7**2 + forceJetsz7**2)


forceTotSlabs = np.sqrt(forceSlabsx**2 + forceSlabsy**2 + forceSlabsz**2)


forceTotJets = forceTotJets1 + forceTotJets2 + forceTotJets3 + forceTotJets4 + forceTotJets5 + forceTotJets6 + forceTotJets7

# Calculate impulse on each row. Dummy data for density height and area of the cylinders because the velocity is not needed for now
# impulseJets = calculateVelocity(forceTotJets , time , rho , H , Area*1.5)[1] # First row is composed by 1.5 cylinders
# impulseSlabs = calculateVelocity(forceTotSlabs , time ,rho , H , Area*4.5)[1] # First row is composed by 4.5 cylinders

# velocityRow1 = np.array(calculateVelocity(forceTotRow1 , time , rho , H , Area*1.5)[0]) # First row is composed by 1.5 cylinders
# velocityRow2 = np.array(calculateVelocity(forceTotRow2 , time ,rho , H , Area*4.5)[0]) # First row is composed by 4.5 cylinders

# averageVelocity = (velocityRow1 + velocityRow2 + velocityRow3 + velocityRow4 + velocityRow5 + velocityRow6 + velocityRow7)/7

# Plot total force
plt.figure(1)
plt.plot(forceTotSlabs, label='Slabs')
plt.plot(forceTotJets, label='Jets')
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Force (N)')
plt.legend()
plt.show()




