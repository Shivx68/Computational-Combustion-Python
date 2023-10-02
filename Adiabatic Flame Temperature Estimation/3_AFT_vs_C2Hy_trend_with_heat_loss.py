"""
Constant pressure combustion of Ethane, Ethene and Ethyne with heat loss = 35% of maximum heat

Stoichiometric combustion of hydrocarbons:

C_xH_y + (x + y/4) (O2 + 3.76 N2) -----> x CO2 + (y/2) H2O + 3.76*(x + y/4) N2
"""
import matplotlib.pyplot as plt

# coefficients data for calculating enthalpy
# reactants at 298.15 K (<1000 K)
c2h6_coeffs_r = [4.29142492E+00, -5.50154270E-03, 5.99438288E-05, -7.08466285E-08, 2.68685771E-11, -1.15222055E+04]
c2h4_coeffs_r = [3.95920148E+00, -7.57052247E-03, 5.70990292E-05, -6.91588753E-08, 2.69884373E-11, 5.08977593E+03]
c2h2_coeffs_r = [8.08681094E-01, 2.33615629E-02, -3.55171815E-05, 2.80152437E-08, -8.50072974E-12, 2.64289807E+04]
o2_coeffs_r = [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03]
n2_coeffs_r = [0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04]

# products at AFT (>1000 K)
co2_coeffs_p = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04]
h2o_coeffs_p = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04]
n2_coeffs_p = [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04]

R = 8.314 # J/mol-K
T_std = 298.15 # K
H_loss = 0.35
nh = [6, 4, 2] # no. of hydrogen atoms [ethane, ethene, ethyne]
T_guess = [0]*len(nh)

# function to compute enthalpy
def h(coeffs,T):
	a1 = coeffs[0]
	a2 = coeffs[1]
	a3 = coeffs[2]
	a4 = coeffs[3]
	a5 = coeffs[4]
	a6 = coeffs[5]

	return (a1 + (a2/2)*T + (a3/3)*pow(T,2) + (a4/4)*pow(T,3) + (a5/5)*pow(T,4) + a6/T)*R*T

for i in range(0,len(nh)):
	def f(T):
		# internal energy balance function
		x = 2	# no. of C atoms
		y = nh[i]

		if y is 6:
			fuel_coeffs_r = c2h6_coeffs_r
		if y is 4:
			fuel_coeffs_r = c2h4_coeffs_r
		if y is 2:
			fuel_coeffs_r = c2h2_coeffs_r

		h_reactants = h(fuel_coeffs_r, T_std) + (x + y/4)*(h(o2_coeffs_r, T_std) + 3.76*h(n2_coeffs_r, T_std))
		h_products = x*h(co2_coeffs_p, T) + (y/2)*h(h2o_coeffs_p, T) + (x + y/4)*3.76*h(n2_coeffs_p, T)
		dh_max = h_reactants - x*h(co2_coeffs_p, T_std) + (y/2)*h(h2o_coeffs_p, T_std) + (x + y/4)*3.76*h(n2_coeffs_p, T_std)

		return h_reactants - h_products + 0.35*dh_max # for constant volume combustion, change in internal energy is 0

	def fprime(T):
		# numerical approximation for first order derivative of f(T) using central differencing		
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

print (T_guess)
label = ['C2H6', 'C2H4', 'C2H2']
plt.plot(label,T_guess)
plt.xlabel('Increase in no. of C-C bonds')
plt.ylabel('Temperature (K)')
plt.title('Adiabatic Flame Temperature variation with no. of C-C bonds in\ncombustion of C2 hydrocarbons at constant pressure with 35% heat loss')
plt.show()