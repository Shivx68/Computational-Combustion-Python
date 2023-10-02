"""
Program to study the effect of preheating temperature on the adiabatic flame temperauture at constant volume
"""

import cantera as ct
import matplotlib.pyplot as plt

min_temp = 298
max_temp = 600
T_air = list(range(min_temp, max_temp+1))
plot_interval = 10

gas = ct.Solution('gri30.cti')

for i in T_air:
	
	air = ct.Quantity(gas)
	air.TPX = i, ct.one_atm, 'O2: 0.21, N2: 0.79'
	air.moles = 9.52

	fuel = ct.Quantity(gas)
	fuel.TPX = 298, ct.one_atm, 'CH4: 1'
	fuel.moles = 1

	mix = air + fuel
	print('T_air = ', i)
	print('Temp. before combustion = ', mix.T)

	mix.equilibrate('UV')
	print('T_ad = ', mix.T)
	print('Enthalpy = ', mix.h)

	if i%plot_interval is 0:
		plt.plot(i, mix.T, '*')
	
plt.title('Effect of Preheating on Adiabatic Flame Temperature')
plt.xlabel('Temperature of preheated air at inlet (K)')
plt.ylabel('Adiabatic flame temperature (K)')
plt.show()