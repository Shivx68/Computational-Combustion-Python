"""
To write a Python code using the Cantera module to determine and plot the 'n' number of most sensitive reactions to Temperature out of the 'N' number of total reactions in the GRI30 mechanism for the auto-ignition of methane.
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# total number of reaction parameters being considered
n_param = 100
S_max = [0]*n_param

# number of most sensitive reactions to be plotted
n_plot = 10
S_plot = [0]*n_plot

# time control
t_end = 2e-3
dt = 5e-6

gas = ct.Solution('gri30.cti')
temp = 1500 # K
pres = ct.one_atm # Pa

gas.TPX = temp, pres, 'CH4:1, O2:2, N2:7.52'
r = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([r])

# including all the reaction parameters with respect to which the Temperature sensitivity is computed
for i in range(n_param):
	r.add_sensitivity_reaction(i)

# setting tolerances
sim.rtol = 1e-6
sim.atol = 1e-15
sim.rtol_sensitivity = 1e-6
sim.atol_sensitivity = 1e-6

# time integration loop
for t in np.arange(0, t_end, dt):	
	sim.advance(t)
	
	# computing sensitivities for each reaction at every time step
	for j in range(n_param):		
		S = sim.sensitivity(1, j)
		
		# checking for maximum sensitivity value at each time step for every reaction
		if np.abs(S) > np.abs(S_max[j]):
			S_max[j] = S

label = []

# gleaning the 'n_plot' most sensitive reactions from 'n-param' total reactions
for k in range(n_plot):
	maximum = np.max(S_max)
	minimum = np.min(S_max)

	# storing the most sensitive reactions inside 'S_plot' in the descending order of sensitivities
	if np.abs(maximum) > np.abs(minimum):
		S_plot[k] = maximum
		reaction_index = S_max.index(maximum)
		label.append(sim.sensitivity_parameter_name(reaction_index))
		S_max.remove(maximum)
	else:
		S_plot[k] = minimum
		reaction_index = S_max.index(minimum)
		label.append(sim.sensitivity_parameter_name(reaction_index))
		S_max.remove(minimum)

# printing results	
for i in range(n_plot):
	print(label[i],'\t', S_plot[i])

# deleting/replacing the redundant characters in each reaction label inside string array 'label'
red_chars = 'IdealGasConstPressureReactor_1: '
rep = ''
for r in range(len(label)):
	label[r] = label[r].replace(red_chars, rep)

# decreasing the 'Y' axis ticklabel size
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.tick_params(axis='y', labelsize=8)

# plotting the desired results
plt.title('Most sensitive reaction parameters to Temperature\nin the auto-ignition of CH4')
plt.barh(label, S_plot)
plt.xlabel('Temperature sensitivity')
plt.tight_layout()
plt.show()