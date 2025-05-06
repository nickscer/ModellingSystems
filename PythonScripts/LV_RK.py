import numpy as np
import matplotlib.pyplot as plt

# Parameters
r1 = 0.61      # idealistic growth rate of red squirrel 
r2 = 0.82      # idealistic growth rate of grey squirrel (higher cuz greys breed faster)
'''(a11 higher cuz native habitat and ability to avoid predators like pine martens)'''
a11 = -2.44 * 10 ** -7   # limitation for red squirrel  i.e. yields carrying capacity along with r1  
a12 =  1.952 * 10 ** -6   # competition towards red squirrel  i.e. how much grey squirrel hurts red squirrel (higher cuz SQPV and size)
a21 =  2.46 * 10 ** -8   # competition towards grey squirrel i.e. how much red squirrel hurts grey squirrel 
a22 = -2.73333333 * 10 ** -7   # limitation for grey squirrel i.e. yields carrying capacity along with r2 
u0 = [10000, 28000] # ICs
# Equilibria at [107758.6207, 2990301.724]


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
T_max = 25
t_values = np.arange(0, T_max + h, h)
u = np.zeros([2, len(t_values)])
print(u)

u[:, 0] = u0
for i in range(1, len(t_values)):
    u_out = rk4(squirrely, h, t_values[i - 1], u[:, i - 1])
    u[:, i] = u_out
print(u)

# Plotting
fig, ax = plt.subplots(2,1)
ax[0].plot(t_values, u[0, :], label='Red Squirrels', color='red')
ax[0].plot(t_values, u[1, :], label='Grey Squirrels', color='grey')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Population Size')
ax[0].set_title('Population Dynamics of Red and Grey Squirrels')
ax[0].legend()

ax[1].plot(u[0, :], u[1, :])
ax[1].set_xlabel('Red Squirrels')
ax[1].set_ylabel('Grey Squirrels')

fig.tight_layout()
plt.show()