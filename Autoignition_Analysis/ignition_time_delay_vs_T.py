"""
Calculating the ignition time delay for auto-ignition of methane at different temperatures
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np


t_end = 10
dt = 1e-4
T_start = 950
T_end = 1450
dT = 1
P = 5*ct.one_atm
t_delay = []

T = np.arange(T_start, T_end+dT, dT)
gas = ct.Solution('gri30.cti')


for i in range(len(T)):
	
	# calculating the ignition delay temperature
	T_ign = T[i] + 400
	gas.TPX = T[i], P, 'CH4:1, O2:2, N2:7.52'
	
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
plt.plot(T, t_delay)
plt.xlabel('Temperature (K)')
plt.ylabel('Ignition time delay (ms)')
plt.title('Ignition time delay variation with Temperature for the auto-ignition of CH4')
plt.show()




