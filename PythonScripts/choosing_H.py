import numpy as np
import matplotlib.pyplot as plt

# Parameters
r1 = 0.61
r2 = 0.82
a11 = -2.44e-7
a12 = -1.952e-7
a21 = -2.46e-8
a22 = -2.73333333e-7

u0 = np.array([3.5e6, 1e4])
t0, T = 1930, 2025
dt = 0.1
t_values = np.arange(t0, T + dt, dt)

# RK4
def rk4(fun, dt, t, u):
    k1 = fun(t, u)
    k2 = fun(t + dt/2, u + dt/2 * k1)
    k3 = fun(t + dt/2, u + dt/2 * k2)
    k4 = fun(t + dt,   u + dt * k3)
    return u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

def squirrely(t, u, H):
    R, G = u
    dR = R * (r1 + a11 * R + a12 * G)
    dG = G * (r2 + a21 * R + a22 * G + H)
    return np.array([dR, dG])

# Heuristic H values to test
H_values = [-0.2, -0.1, -0.05, -0.03, -0.01]

fig, (axR, axG) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

R_target, G_target = 140000, 2500000
axR.axhline(R_target, linestyle='--', label='Red target')
axG.axhline(G_target, linestyle='--', label='Grey target')

for H in H_values:
    u = np.zeros((2, len(t_values)))
    u[:, 0] = u0
    for i in range(1, len(t_values)):
        u[:, i] = rk4(lambda t, y: squirrely(t, y, H), dt, t_values[i-1], u[:, i-1])
    axR.plot(t_values, u[0], label=f'H={H}')
    axG.plot(t_values, u[1], label=f'H={H}')

# Plotting
axR.set_ylabel('Red population')
axG.set_ylabel('Grey population')
axG.set_xlabel('Year')
axR.legend()
axG.legend()
plt.tight_layout()
plt.show()
