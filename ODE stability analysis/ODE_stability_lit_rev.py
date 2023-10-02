"""
Stability analysis for the ODE
dy/dt = -1000*y + 3000 - 2000*e^-t

The analytical solution is
y = 3 - 0.998*e^(-1000*t) - 2.002*e^-t

The numerical explicit solution is
y_(n+1) = y_n + h*(-1000*y_n + 3000 - 2000*e^-t)
"""

import numpy as np
import matplotlib.pyplot as plt

t_end = 0.5	# simulation time
h = 0.002	# step size
t_array = np.arange(h, t_end, h)
time = np.arange(0, t_end, h)
y0 = 0	# initial value
y_numerical = [y0]
y_analytical = [y0]
y_old = y0

# loop for solution advance
for t in t_array:
	
	# numerical solution using explicit Euler or forward difference
	y_new = y_old + h*(-1000*y_old + 3000 - 2000*np.exp(-t))
	y_numerical.append(y_new)
	y_old = y_new

	# analytical solution
	y_analytical.append(3 - 0.998*np.exp(-1000*t) - 2.002*np.exp(-t))

# plot
fig1 = plt.figure(1)
plt.plot(time, y_numerical, alpha = 0.7, label = 'Numerical Sol')
plt.plot(time, y_analytical, alpha = 0.7, label = 'Analytical Sol')
plt.legend(loc = 'best')
plt.xlabel('time')
plt.ylabel('y = f(t)')
fig1.suptitle('Stability analysis for a simple ODE for time step h = {}'.format(h))
plt.show()