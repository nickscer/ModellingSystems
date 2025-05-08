import numpy as np
import matplotlib.pyplot as plt

# Parameters
r1 = 0.61   # idealistic growth rate of red squirrel 
r2 = 0.82   # idealistic growth rate of grey squirrel (higher cuz greys breed faster)
'''(a11 higher cuz native habitat and ability to avoid predators like pine martens)'''
a11 = -2.44 * 10 ** -7    # limitation for red squirrel  i.e. yields carrying capacity along with r1  
a12 = -1.952 * 10 ** -7   # competition towards red squirrel  i.e. how much grey squirrel hurts red squirrel (higher cuz SQPV and size)
a21 = -2.46 * 10 ** -8    # competition towards grey squirrel i.e. how much red squirrel hurts grey squirrel 
a22 = -2.73333333 * 10 ** -7   # limitation for grey squirrel i.e. yields carrying capacity along with r2 
u0 = [3500000, 1000] # ICs
# Equilibria at [107758, 2990517]


# Functions setup
'''xdot = x(r1 + a11*x + a12*y)'''
'''ydot = y(r2 + a21*x + a22*y)'''

def rk4(fun, dt, t0, u0):
    """Runge-Kutta 4th order method
    fun: function that returns the derivative of u at time t
    """
    k1 = fun(t0, u0)
    k2 = fun(t0 + dt / 2, u0 + dt / 2 * k1)
    k3 = fun(t0 + dt / 2, u0 + dt / 2 * k2)
    k4 = fun(t0 + dt, u0 + dt * k3)
    u_out = u0 + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return u_out

def squirrely(t0, u):
    """Function that returns the derivative of u at time t
    u is array s.t. [x, y] where x is reds population and y is greys population
    """
    udot = [u[0] * (r1 + a11 * u[0] + a12 * u[1]), 
         u[1] * (r2 + a21 * u[0] + a22 * u[1])]
    return np.array(udot)

# Computing values 
h = 0.001
T_max = 115 # 115 years
t_values = np.arange(0+1930, T_max + 1930 + h, h)
u = np.zeros([2, len(t_values)])
print(u)

u[:, 0] = u0
for i in range(1, len(t_values)):
    u_out = rk4(squirrely, h, t_values[i - 1], u[:, i - 1])
    u[:, i] = u_out
print(u)

# Plotting the population dynamics (ax[0]) in a separate figure
fig1, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(t_values, u[0, :], label='Red Squirrels', color='red')
ax1.plot(t_values, u[1, :], label='Grey Squirrels', color='grey')
ax1.set_xlabel('Time in years starting from 1930')
ax1.set_ylabel('Population Size')
ax1.set_title('Population Dynamics of Red and Grey Squirrels')
ax1.axhline(107758, color='red', linestyle='--', label='Red Equilibrium')
ax1.axhline(2990301, color='grey', linestyle='--', label='Grey Equilibrium')
ax1.legend()
fig1.tight_layout()
plt.show()

# Define a grid around the stable saddle point
x_vals = np.linspace(107758 - 50000, 107758 + 50000, 20)  # Red squirrel range
y_vals = np.linspace(2990301 - 500000, 2990301 + 500000, 20)  # Grey squirrel range
X, Y = np.meshgrid(x_vals, y_vals)

# Compute the vector field
U = X * (r1 + a11 * X + a12 * Y)  # dx/dt
V = Y * (r2 + a21 * X + a22 * Y)  # dy/dt

# Add a new subplot for the phase portrait with vectors using quiver
# Define a grid covering the entire trajectory range
x_min, x_max = 0, 4e6  # Adjust based on trajectory data
y_min, y_max = -0.1e6, 3.5e6
x_vals = np.linspace(x_min, x_max, 25)
y_vals = np.linspace(y_min, y_max, 25)
X, Y = np.meshgrid(x_vals, y_vals)

# Compute the vector field
U = X * (r1 + a11 * X + a12 * Y)
V = Y * (r2 + a21 * X + a22 * Y)

# Plotting phase space(s)
fig2, (ax2, ax3) = plt.subplots(1, 2, figsize=(16, 8))

# Phase space
ax2.plot(u[0, :], u[1, :], label='Trajectory', color='black', alpha=0.7)
ax2.set_xlabel('Red Squirrels')
ax2.set_ylabel('Grey Squirrels')
ax2.set_title('Phase Space of Red and Grey Squirrels')
ax2.plot(107758, 2990301, 'go', markersize=10, label='Stable sink point $(107758, 2990301)$')
ax2.legend()
ax2.set_xlim(x_min, x_max)  # Synchronize x-axis
ax2.set_ylim(y_min, y_max)  # Synchronize y-axis
ax2.set_aspect('equal', adjustable='box')  # Equal aspect ratio

# Phase space with vector field
ax3.plot(u[0, :], u[1, :], label='Trajectory', color='black', alpha=0.7)
ax3.quiver(X, Y, U, V, color='blue', angles='xy', 
           scale_units='xy', scale=1.7, width=0.002, headwidth=3)
ax3.set_xlabel('Red Squirrels')
ax3.set_ylabel('Grey Squirrels')
ax3.set_title('Phase Portrait with Vector Field')
ax3.plot(107758, 2990301, 'go', markersize=10, label='Stable sink point')
ax3.set_xlim(x_min, x_max)  # Synchronize x-axis
ax3.set_ylim(y_min, y_max)  # Synchronize y-axis
ax3.set_aspect('equal', adjustable='box')  # Equal aspect ratio
ax3.legend()

fig2.tight_layout()
plt.show()