"""
Program to study the effect of preheating on the effeciency of combustion of methane in air
"""

import cantera as ct
import matplotlib.pyplot as plt

phi = 1
min_temp = 298
max_temp = 600
T_air = list(range(min_temp, max_temp+1))
plot_interval = 10
LHV = 50e6 # J/kg

gas = ct.Solution('gri30.cti')

# defining fuel quantity
fuel = ct.Quantity(gas)
fuel.TPX = 298, ct.one_atm, 'CH4: 1'
fuel.moles = 1
energy_fuel = fuel.h * fuel.mass / 1000

# defining exhaust gas quantity
flue = ct.Quantity(gas)
flue.TPX = 1700, ct.one_atm, 'CO2:1, H2O:2, N2:7.52'
flue.moles = 1 + 2 + 7.52
energy_flue = flue.h * flue.mass / 1000

# loop to iterate air inlet temperature
for i in T_air:

    # defining air quantity
    air = ct.Quantity(gas)
    air.TPX = i, ct.one_atm, 'O2: 0.21, N2: 0.79'
    air.moles = 9.52 / phi
    m_air = 28.84 * air.moles / 1000
    energy_air = air.h * air.mass / 1000

    # energy balance
    Q_out = (energy_air  + energy_fuel) - energy_flue
    Q_out_sp = Q_out / (fuel.mass * 0.001)
    efficiency = Q_out_sp / LHV
    print('T_air = ', i, 'K\tEfficiency = ', efficiency, '%')
    
    if i % plot_interval is 0:
        plt.plot(i, efficiency,'o') 

plt.title('Effect of preheating on the efficiency of combustion of CH4 in air')
plt.xlabel('Temperature of preheated air at inlet (K)')
plt.ylabel('Efficiency')
plt.show()
