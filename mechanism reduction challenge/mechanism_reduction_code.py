"""
GRI30 mechanism reduction based on temperature sensitivity

"""

import os
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Extracting input date from 'input_file' using user-defined class 'FileReader'
class FileReader:

	def __init__(self, filename):
		self.filename = filename
		self.input_dict = {}

	def validate_file(self):
		return os.path.isfile(self.filename)

	def read(self):
		if self.validate_file():
			for line in open(self.filename, 'r'):
				self.input_dict[line.split()[0]] = line.split()[1:]
				
		else:
			print('File name does not exist')

a = FileReader('input.txt')
a.read()
T_val = [float(x) for x in a.input_dict['Temperature']]
P_val = [float(x) for x in a.input_dict['Pressure']]
phi_val = [float(x) for x in a.input_dict['Phi']]
for x in a.input_dict['End_time']: t_end = float(x)
for x in a.input_dict['Delta_t']: dt = float(x)
for x in a.input_dict['Min_size']: min_size = int(x)
for x in a.input_dict['Tol_Tmax']: tol_T_max = float(x)
for x in a.input_dict['Tol_igd']: tol_ign_delay = float(x)
for x in a.input_dict['Fuel']: fuel = x

#################################################

# Ordering the reaction indices in decreasing order of sensitivity to temperature for all combination of state values

S = []
S_max0 = []
R_index = []
const1 = True
n = 0

# air mixture 
air = 'O2:1, N2:3.76'

# creating gas object
gas = ct.Solution('gri30.cti')
S_max = [0]*gas.n_reactions

# state values combination loop
for T in T_val:

	for P in P_val:

		for phi in phi_val:
			
			gas.TP = T, P
			gas.set_equivalence_ratio(phi, fuel, air)

			# creating a reactor and reactor net objcets
			r = ct.IdealGasReactor(gas)
			sim = ct.ReactorNet([r])

			# including all reactions from gri mech to perform sensitivity analysis
			for i in range(0, gas.n_reactions):
				r.add_sensitivity_reaction(i)

			# setting absolute and relative tolerances
			sim.rtol = 1e-6
			sim.atol = 1e-15
			sim.rtol_sensitivity = 1e-6
			sim.atol_sensitivity = 1e-6
		

			# solving the reaction rate odes and sensitivity odes
			for t in np.arange(0, t_end+dt, dt):

				sim.advance(t)

				# finding maximum sensitivities for all reactions
				for i in range(0, gas.n_reactions):
					
					S = sim.sensitivity(2, i)

					if np.abs(S) > np.abs(S_max[i]):
						S_max[i] = np.abs(S)

			# zipping the max sensitivities with respective reaction indices
			for i in range(0, gas.n_reactions):
				S_max0.append([])
				S_max0[i].append(np.abs(S_max[i]))
				S_max0[i].append(i)

			# sorting S_max0 array in decreasing order of sensitivity
			S_sorted = sorted(S_max0, key = lambda x: x[0], reverse = True)
			
			# appending the sorted order of reaction indices to R_index
			for i in range(0, gas.n_reactions):

				if const1 is True:
					R_index.append([])

				R_index[i].append(S_sorted[i][1])
			
			# resetting values
			const1 = False
			S_max0.clear()
			n = n + 1

#################################################

# flattening out R_index to create a single list of reactions from all state value combinations in the order of decreasing sensitivities

R_set = set()
R_ordered = []

for i in range(0, gas.n_reactions):

	for j in range(0, n):

		if R_index[i][j] not in R_set:
			R_ordered.append(R_index[i][j])
			R_set.add(R_index[i][j])

print(R_ordered,'\n')
#################################################

# Determining the size of the reduced mechanism

spec = ct.Species.listFromFile('gri30.cti')
rxns = []
ign_delay_arr = []
T_max_arr = []
mech_size = []

for T in T_val:
	for P in P_val:
		for phi in phi_val:
			
			# computing reference ign_delay and T_max for each state
			T_ign = T + 400	# K

			gas.TP = T, P
			gas.set_equivalence_ratio(phi, fuel, air)

			r = ct.IdealGasReactor(gas)
			sim = ct.ReactorNet([r])
			
			const2 = True
			T_max_ref = 0

			for t in np.arange(0, t_end+dt, dt):

				sim.advance(t)

				if (gas.T > T_ign and const2 is True):
					ign_delay_ref = t
					const2 = False
				if (gas.T > T_max_ref):
					T_max_ref = gas.T

			# performing mechanism reduction for each state
			err_T_max = 100
			err_ign_delay = 100
			i = 0
			const4 = True
			
			# appending reactions one by one and comparing T_max and ign_delay results
			for i in range(0, gas.n_reactions):

				rxns.append(gas.reaction(R_ordered[i]))
				gas1 = ct.Solution(thermo = 'IdealGas', kinetics = 'GasKinetics', species = spec, reactions = rxns)
				gas1.TP = T, P
				gas1.set_equivalence_ratio(phi, fuel, air)
				r1 = ct.IdealGasReactor(gas1)
				sim1 = ct.ReactorNet([r1])
				const3 = True
				T_max = 0

				for t in np.arange(0, t+dt, dt):

					sim1.advance(t)

					if (gas1.T > T_ign and const3 is True):
						ign_delay = t
						const3 = False

					if (gas1.T > T_max):
						T_max = gas1.T

				# incrementing reaction index
				i = i + 1

				# catching exception until the no. of reactions increase to suitable number
				try :
					
					# storing values for the plot
					if (i >= min_size) :
						ign_delay_arr.append(ign_delay)
						T_max_arr.append(T_max)
						mech_size.append(i)

					# computing accuracy
					err_T_max = np.abs(((T_max_ref - T_max)/T_max_ref)*100) # % error in maximum temperature
					err_ign_delay = np.abs(((ign_delay_ref - ign_delay)/ign_delay_ref)*100)	# % error in ignition delay

					# printing mechanism size and other results when the tolerance criteria is satisfied
					if (err_T_max < tol_T_max and err_ign_delay < tol_ign_delay and const4 is True):
						
						print('T = %g, P = %g, phi = %g, T_max_ref = %g, ign_delay_ref = %g, T_max = %g, ign_delay = %g, T_max_error = %g, Ign_delay_error = %g, red_mech_size = %g' %(T, P, phi, T_max_ref, ign_delay_ref, T_max, ign_delay, err_T_max, err_ign_delay, i))
						const4 = False

				except NameError:

					pass

			
			# plotting results
			plt.subplot(2,1,1)
			plt.title('State Values: T = %g K, P = %g bar, phi = %g' %(T, P/100000, phi))
			plt.plot(mech_size, ign_delay_arr, color = 'green')
			plt.axhline(y = ign_delay_ref, color = 'red', linestyle = 'dotted', label = 'Reference Ignition Delay')
			plt.xlabel('No. of reactions most sensitive to Temperature')
			plt.ylabel('Ignition Delay (s)')
			plt.legend(loc = 'best')

			plt.subplot(2,1,2)
			plt.plot(mech_size, T_max_arr, color = 'blue')
			plt.axhline(y = T_max_ref, color = 'red', linestyle = 'dotted', label = 'Reference Max Temperature')
			plt.xlabel('No. of reactions most sensitive to Temperature')
			plt.ylabel('Maximum Temperature (K)')
			plt.legend(loc = 'best')

			plt.tight_layout()
			plt.show()

			# resetting arrays for next state values combination
			rxns.clear()
			ign_delay_arr.clear()
			T_max_arr.clear()
			mech_size.clear()