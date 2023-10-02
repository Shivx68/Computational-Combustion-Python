"""
Constant pressure combustion of methane  with heat loss = 35% of maximum heat

Stoichiometric combustion of methane:
CH4 + 2 (O2 + 3.76 N2) -----> CO2 + 2 H2O + 7.52 N2
"""
import matplotlib.pyplot as plt

# coefficients data for calculating enthalpy
# reactants at 298.15 K (<1000 K)
ch4_coeffs_r = [5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04]
o2_coeffs_r = [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03]
n2_coeffs_r = [0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04]

# products at AFT (>1000 K)
co2_coeffs_p = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04]
h2o_coeffs_p = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04]
n2_coeffs_p = [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04]

R = 8.314 # J/mol-K
T_std = 298.15 # K
H_loss = 0.35

# function to compute enthalpy
def h(coeffs,T):
	a1 = coeffs[0]
	a2 = coeffs[1]
	a3 = coeffs[2]
	a4 = coeffs[3]
	a5 = coeffs[4]
	a6 = coeffs[5]

	return (a1 + (a2/2)*T + (a3/3)*pow(T,2) + (a4/4)*pow(T,3) + (a5/5)*pow(T,4) + a6/T)*R*T

def f(T):
		# internal energy balance function		
		
		h_reactants = h(ch4_coeffs_r, T_std) + 2*(h(o2_coeffs_r, T_std) + 3.76*h(n2_coeffs_r, T_std))
		h_products = h(co2_coeffs_p, T) + 2*h(h2o_coeffs_p, T) + 7.52*h(n2_coeffs_p, T)
		dh_max = h_reactants - h(co2_coeffs_p, T_std) + 2*h(h2o_coeffs_p, T_std) + 7.52*h(n2_coeffs_p, T_std)

		return h_reactants - h_products + 0.35*dh_max # for constant volume combustion, change in internal energy is 0

def fprime(T):
	# numerical approximation for first order derivative of f(T) using central differencing		
	dT = 1e-6
	return (f(T+dT) - f(T-dT))/(2*dT)

# Newton-Raphson loop
T_guess = 1500
tol = 1e-3
alpha = 0.5
count = 1

while abs(f(T_guess)) > tol:
	T_guess = T_guess - alpha*(f(T_guess)/fprime(T_guess))
	count = count + 1

print (T_guess)