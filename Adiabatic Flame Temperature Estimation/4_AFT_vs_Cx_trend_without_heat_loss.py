"""
Constant pressure combustion of various alkanes with no heat loss

General stoichiometric reaction for alkane combustion

C(n)H(2n+2) + (3n+1)/2 (O2 + 3.76 N2) -----> n CO2 + (n+1) H20 + 3.76*(3n+1)/2 N2
"""

import matplotlib.pyplot as plt
import cantera as ct
gas = ct.Solution('gri30.cti')

# coefficients data for calculating enthalpy
# reactants at 298.15 K (<1000 K)
ch4_coeffs_r = [5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04]
c2h6_coeffs_r = [4.29142492E+00, -5.50154270E-03, 5.99438288E-05, -7.08466285E-08, 2.68685771E-11, -1.15222055E+04]
c3h8_coeffs_r = [0.93355381E+00, 0.26424579E-01, 0.61059727E-05, -0.21977499E-07, 0.95149253E-11, -0.13958520E+05]
o2_coeffs_r = [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03]
n2_coeffs_r = [0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04]

# products at AFT (>1000 K)
co2_coeffs_p = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04]
h2o_coeffs_p = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04]
n2_coeffs_p = [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04]

R = 8.314 # J/mol-K
T_std = 298.15 # K
T_guess = [0]*3
T_ad = [0]*3
# function to compute enthalpy
def h(coeffs,T):
	a1 = coeffs[0]
	a2 = coeffs[1]
	a3 = coeffs[2]
	a4 = coeffs[3]
	a5 = coeffs[4]
	a6 = coeffs[5]

	return (a1 + (a2/2)*T + (a3/3)*pow(T,2) + (a4/4)*pow(T,3) + (a5/5)*pow(T,4) + a6/T)*R*T

for n in range(1,4):
	def f(T):
		# internal energy balance function		
		if n is 1:
			fuel_coeffs_r = ch4_coeffs_r
		if n is 2:
			fuel_coeffs_r = c2h6_coeffs_r
		if n is 3:
			fuel_coeffs_r = c3h8_coeffs_r

		h_reactants = h(fuel_coeffs_r, T_std) + (3*n+1)/2 *(h(o2_coeffs_r, T_std) + 3.76*h(n2_coeffs_r, T_std))
		h_products = n*h(co2_coeffs_p, T) + (n+1)*h(h2o_coeffs_p, T) + 3.76*(3*n+1)/2 *h(n2_coeffs_p, T)

		return h_reactants - h_products

	def fprime(T):
		# numerical approximation for first order derivative of f(T) using central differencing		
		dT = 1e-6
		return (f(T+dT) - f(T-dT))/(2*dT)

	# Newton-Raphson loop 
	T_guess[n-1] = 1500
	tol = 1e-3
	alpha = 0.5
	count = 1

	while abs(f(T_guess[n-1])) > tol:
		T_guess[n-1] = T_guess[n-1] - alpha*(f(T_guess[n-1])/fprime(T_guess[n-1]))
		count = count + 1

# Computing AFT for alkanes using Cantera
gas = ct.Solution('gri30.cti')

n=1
n_ch4 = 1
n_o2 = (3*n+1)/2
n_n2 = 3.76*(3*n+1)/2
n_tot = n_ch4 + n_o2 + n_n2
gas.TPX = 298.15, 101325, {'CH4': n_ch4/n_tot, 'O2': n_o2/n_tot, 'N2': n_n2/n_tot}
gas.equilibrate('HP')
T_ad[0] = gas.T

n=2
n_c2h6 = 1
n_o2 = (3*n+1)/2
n_n2 = 3.76*(3*n+1)/2
n_tot = n_c2h6 + n_o2 + n_n2
gas.TPX = 298.15, 101325, {'C2H6': n_c2h6/n_tot, 'O2': n_o2/n_tot, 'N2': n_n2/n_tot}
gas.equilibrate('HP')
T_ad[1] = gas.T

n=3
n_c3h8 = 1
n_o2 = (3*n+1)/2
n_n2 = 3.76*(3*n+1)/2
n_tot = n_c3h8 + n_o2 + n_n2
gas.TPX = 298.15, 101325, {'C3H8': n_c3h8/n_tot, 'O2': n_o2/n_tot, 'N2': n_n2/n_tot}
gas.equilibrate('HP')
T_ad[2] = gas.T

print(T_guess)
print(T_ad)
label = ['CH4', 'C2H6', 'C3H8']
plt.plot(label, T_guess, label = 'Without Cantera')
plt.plot(label, T_ad, label = 'With Cantera')
plt.xlabel('Increase in no. of C-atoms')
plt.ylabel('Temperature (K)')
plt.legend(loc = 'upper left')
plt.title('Adiabatic Flame Temperature variation with no. of C- atoms in\ncombustion of alkanes at constant pressure with no heat loss')
plt.show()