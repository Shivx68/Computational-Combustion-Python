"""
flame speed analysis for the combustion of hydrogen
"""

import cantera as ct
import matplotlib.pyplot as plt

# initial conditions
temp = 300
pres = ct.one_atm
reactants = 'H2:2, O2:1, N2:3.76'

# defining the gas object
gas =  ct.Solution('h2_mech.cti')
gas.TPX = temp, pres, reactants

# defining the domain width for 1d flame peopagation
width = 0.03

# creating a freeflame object
f = ct.FreeFlame(gas, width = width)

# setting the refine criteria for automated geid refinement
f.set_refine_criteria(ratio = 3, slope = 0.07, curve = 0.14)
f.solve()
conc = f.concentrations
index = gas.species_index('CO2')


# plotting the results
fig = plt.figure(1)
plt.subplot(2,1,1)
plt.plot(f.grid, f.T)
plt.ylabel('Temperature (K)')
plt.xlabel('Domain length (m)')
plt.title('Flame Speed analysis for combustion of Hydrogen')
plt.subplot(2,1,2)
plt.plot(f.grid, f.u)
plt.ylabel('Velocity (m/s)')
plt.xlabel('Domain length (m)')
fig.tight_layout()

fig2 = plt.figure(2)
plt.plot(f.grid, conc[index,:])
plt.xlabel('Domain length (m)')
plt.ylabel('Concentration of CO2 (kmol/m^3)')
fig2.suptitle('Variation in concentration of CO2 across the domain')
plt.show()