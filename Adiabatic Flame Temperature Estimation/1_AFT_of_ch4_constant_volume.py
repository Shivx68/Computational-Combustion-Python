"""
Computing the Adiabatic Flame Temperature for Constant-Volume combustion of Methane and Studying the effect of equivalent ratio on AFT and comparing the results with those obtained from Cantera

Combustion reaction for equivalence ratio (phi) <=1 :
CH4 + 2/phi (O2 + 3.76 N2) -----> CO2 + 2 H2O + (2/phi - 2) O2 + 7.52/phi N2

Combustion reaction for phi > 1 :
CH4 + O2 + 3.76 N2 -----> 0.5 CO2 + 0.5 H2O + 0.5 CO + 1.5 H2 + 3.76 N2

"""

import matplotlib.pyplot as plt
import cantera as ct

# coefficients data for calculating enthalpy
# reactants at 298.15 K (<1000 K)
ch4_coeffs_r = [5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04]
o2_coeffs_r = [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03]
n2_coeffs_r = [0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04]

# products at AFT (>1000 K)
co2_coeffs_p = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04]
h2o_coeffs_p = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04]
n2_coeffs_p = [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04]
o2_coeffs_p = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03]
co_coeffs_p = [2.71518561E+00, 2.06252743E-03, -9.98825771E-07, 2.30053008E-10, -2.03647716E-14, -1.41518724E+04]
h2_coeffs_p = [ 3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02]


R = 8.314 # J/mol-K
T_std = 298.15
phi = [0.2, 0.4, 0.6, 0.8, 1, 2] # equivalemce ratio
T_guess = [0]*len(phi)
T_ad = [0]*len(phi)

# function to compute enthalpy
def h(coeffs,T):
	a1 = coeffs[0]
	a2 = coeffs[1]
	a3 = coeffs[2]
	a4 = coeffs[3]
	a5 = coeffs[4]
	a6 = coeffs[5]

	return (a1 + (a2/2)*T + (a3/3)*pow(T,2) + (a4/4)*pow(T,3) + (a5/5)*pow(T,4) + a6/T)*R*T

# loop to calculate AFT for different equivalence ratios
for i in range(0,len(phi)):
	def f(T):
		# internal energy balance function		
		
		if phi[i] <= 1:
			h_reactants = h(ch4_coeffs_r, T_std) + (2/phi[i])*h(o2_coeffs_r, T_std) + (7.52/phi[i])*h(n2_coeffs_r, T_std)
			h_products = h(co2_coeffs_p, T) + 2*h(h2o_coeffs_p, T) + (7.52/phi[i])*h(n2_coeffs_p, T) + (2/phi[i] -2)*h(o2_coeffs_p, T)
			dpv = (1 + 2/phi[i] + 7.52/phi[i])*R*(T_std - T) # Const. vol. work, (P_reac - P_prod)*V = R*(n_reac*T_reac - n_prod*T_prod)
		
		else:
			h_reactants = h(ch4_coeffs_r, T_std) + h(o2_coeffs_r, T_std) + 3.76*h(n2_coeffs_r, T_std)
			h_products = 0.5*h(co2_coeffs_p, T) + 0.5*h(h2o_coeffs_p, T) + 3.76*h(n2_coeffs_p, T) + 0.5*h(co_coeffs_p, T) + 1.5*h(h2_coeffs_p, T)
			dpv = R*(5.76*T_std - 6.76*T)

		return h_reactants - h_products - dpv # for constant volume combustion, change in internal energy is 0

	def fprime(T):
		
		dT = 1e-6
		return (f(T+dT) - f(T-dT))/(2*dT)

	# Newton-Raphson loop 
	T_guess[i] = 1500
	tol = 1e-3
	alpha = 0.5
	count = 1

	while abs(f(T_guess[i])) > tol:
		T_guess[i] = T_guess[i] - alpha*(f(T_guess[i])/fprime(T_guess[i]))
		count = count + 1

# Computing Adiabatic Flame Temperature using Cantera
gas = ct.Solution('gri30.xml')

for i in range(0,len(phi)):

	n_ch4 = 1
	n_o2 = 2/phi[i]
	n_n2 = 7.52/phi[i]
	n_tot = n_ch4 + n_o2 + n_n2

	gas.TPX = 298.15, 101325, {'CH4': n_ch4/n_tot, 'O2': n_o2/n_tot, 'N2': n_n2/n_tot}
	gas.equilibrate('UV')
	T_ad[i] = gas.T

print(T_guess)
print(T_ad)

# comparing data with plots
plt.plot(phi,T_guess, color = 'blue', label = 'Without using Cantera')
plt.plot(phi,T_ad, color = 'red', label = 'Using Cantera')
plt.xlabel('Equivalence ratio ($\phi$)')
plt.ylabel('Temperature (K)')
plt.legend(loc = 'upper right')
plt.title('Adiabatic Flame Temperature variation with equivalent ratio in\ncombustion of methane at constant volume')
plt.show()