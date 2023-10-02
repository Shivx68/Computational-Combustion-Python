"""
Calculating the ignition time delay for auto-ignition of methane for different pressures
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# initialization
t_end = 10
dt = 1e-3
T = 1250
T_ign = T + 400
P_start = 1
P_end = 5
dp = 0.01
t_delay = []

P = np.arange(P_start, P_end+dp, dp)
gas = ct.Solution('gri30.cti')

# loop for iterating pressure
for i in range(len(P)):
	
	gas.TPX = T, P[i]*ct.one_atm, 'CH4:1, O2:2, N2:7.52'
	
	# creating a reactor and a reactor nework
	r = ct.IdealGasReactor(gas)
	sim = ct.ReactorNet([r])
	
	# initializing time and T_stop to 0 at start of every loop
	time = 0.0
	T_stop = 0
	
	# convergence loop
	while T_stop < T_ign:

		time = time + dt # in s
		sim.advance(time)
		T_stop = r.T
	
	# storing ignition delay times
	t_delay.append(time*1e3) # in ms

# plotting the trend
plt.plot(P, t_delay)
plt.xlabel('Pressure (atm)')
plt.ylabel('Ignition time delay (ms)')
plt.title('Ignition time delay variation with Pressure for the auto-ignition of CH4')
plt.show()





