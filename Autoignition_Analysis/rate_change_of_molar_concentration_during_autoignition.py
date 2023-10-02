"""
Rate of change of molar concentrations of H2O, O2 and OH at 500 K and 1000 K
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

T_array = [500, 1000]
P = 5*ct.one_atm
t_end = 10
dt = 1e-3
time = np.arange(0, t_end+dt, dt)

gas = ct.Solution('gri30.cti')

for T in T_array:
	gas.TPX = T, P, ('CH4:1, O2:2, N2:7.52')

	r = ct.IdealGasReactor(gas)
	sim = ct.ReactorNet([r])
	sol = ct.SolutionArray(gas, extra = ['time_ms'])

	for i in time:

		sim.advance(i)
		sol.append(r.thermo.state, time_ms = i*1e3)

	plt.plot(time, 10.52*sol.X[:,gas.species_index('H2O')], alpha = 0.5, label = '[H2O]')
	plt.plot(time, 10.52*sol.X[:,gas.species_index('O2')], alpha = 0.5, label = '[O2]')
	plt.plot(time, 10.52*sol.X[:,gas.species_index('OH')], alpha = 0.5, label = '[OH]')
	plt.xlabel('Time (s)')
	plt.ylabel('Molar concentration (moles)')
	title_text = 'Rate of change of molar concentrations of H2O, O2 and OH\nfor auto-ignition of CH4 at {0} K'.format(T)
	plt.title(title_text)
	plt.legend(loc = 'center right')
	plt.show()
