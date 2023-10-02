"""
flame speed analysis for combustion of methane
"""

import cantera as ct
import matplotlib.pyplot as plt

# initial conditions
pres = ct.one_atm
temp = 300.0
reactants = 'CH4:1, O2:2, N2:7.52'

# domain length for 1d flame propagation
width = 0.03

# initializing and defining gas object
gas = ct.Solution('gri30.cti', 'gri30_mix')
gas.TPX = temp, pres, reactants

# creating a freeflame object and passing gas object and width
f = ct.FreeFlame(gas, width = width)

# setting the refine criteria for automated grid refinement
f.set_refine_criteria(ratio = 3, slope = 0.07, curve = 0.14)
f.solve()
conc = f.concentrations
index = gas.species_index('CO2')


# potting the results
fig1 = plt.figure(1)
plt.subplot(2,1,1)
plt.plot(f.grid, f.T)
plt.ylabel('Temperature (K)')
plt.xlabel('Domain length (m)')
plt.title('Flame Speed analysis for combustion of Methane')
plt.subplot(2,1,2)
plt.plot(f.grid, f.u)
plt.ylabel('Velocity (m/s)')
plt.xlabel('Domain length (m)')
fig1.tight_layout()

fig2 = plt.figure(2)
plt.plot(f.grid, conc[index,:])
plt.xlabel('Domain length (m)')
plt.ylabel('Concentration of CO2 (kmol/m^3)')
fig2.suptitle('Variation in concentration of CO2 across the domain')
plt.show()