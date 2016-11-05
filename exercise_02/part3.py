import numpy as np
from plot import frame
from matplotlib import pyplot as plt

# Presets
l1 = 2
l2 = 1
m1 = 0.5
m2 = 1
g = 1
dt = 0.05

iteration_steps = 2000

# Runge-Kutta integrators
def rk2(f, y, t, dt):
    k1 = np.array([ode(t, *y) for ode in f])
    updated_y = y + dt * k1
    k2 = np.array([ode(t + dt, *updated_y) for ode in f])
    return y + 0.5 * dt * (k1 + k2)
def rk4(f, y, t, dt):
    k1 = np.array([ode(t, *y) for ode in f])
    updated_y = y + 0.5 * dt * k1
    k2 = np.array([ode(t + 0.5 * dt, *updated_y) for ode in f])
    updated_y = y + 0.5 * dt * k2
    k3 = np.array([ode(t + 0.5*dt, *updated_y) for ode in f])
    updated_y = y + dt * k3
    k4 = np.array([ode(t + dt, *updated_y) for ode in f])
    return y + dt * (k1/6 + k2/3 + k3/3 + k4/6)

# total energy and relative energy error of y
def total_energy(t, phi1, phi2, q1, q2):
    dphi1 = f_phi1(t, phi1, phi2, q1, q2)
    dphi2 = f_phi2(t, phi1, phi2, q1, q2)
    pot = -(m1 + m2) * g * l1 * np.cos(phi1) - m2 * g * l2 * np.cos(phi2)
    kin = 0.5 * m1 * l1**2 * dphi1**2 + 0.5 * m2 * (l1**2 * dphi1**2 + l2**2 * dphi2**2 + 2*l1*l2*dphi1*dphi2*np.cos(phi1 - phi2))
    return pot + kin
def relative_energy_error(t, all_y):
    return total_energy(t, *[all_y[:, i] for i in range(4)])/total_energy(t, *[all_y[0, i] for i in range(4)]) - 1


# ode's to solve
def f_phi1(t, phi1, phi2, q1, q2):
    return (l2*q1 - l1*q2*np.cos(phi1 - phi2))/(l1**2 * l2 * (m1 + m2 * np.sin(phi1 - phi2)**2))

def f_phi2(t, phi1, phi2, q1, q2):
    return (l1*(m1 + m2)*q2 - l2*m2*q1*np.cos(phi1 - phi2))/(l1 * l2**2 * m2 * (m1 + m2 * np.sin(phi1 - phi2)**2))

def f_q1(t, phi1, phi2, q1, q2):
    return -m2*l1*l2*f_phi1(t, phi1, phi2, q1, q2)*f_phi2(t, phi1, phi2, q1, q2)*np.sin(phi1 - phi2) - (m1 + m2)*g*l1*np.sin(phi1)

def f_q2(t, phi1, phi2, q1, q2):
    return m2*l1*l2*f_phi1(t, phi1, phi2, q1, q2)*f_phi2(t, phi1, phi2, q1, q2)*np.sin(phi1 - phi2) - m2*g*l2*np.sin(phi2)


# start value
y = np.array([np.radians(50), np.radians(-120), 0, 0])

# set f vector
f = [f_phi1, f_phi2, f_q1, f_q2]

t = 0

# save values
all_y = np.array([y])
times = np.array([0])

# integration loop
for i in range(iteration_steps):

    # to use rk2, change function to rk2()
    t += dt
    y = rk4(f, y, t, dt)

    all_y = np.append(all_y, [y], axis=0)
    times = np.append(times, t)

    # if you want to plot, comment this out
    #frame(y[0], y[1], all_y[:, 0], all_y[:, 1], t, i)



### Plot relative energy error
fig = plt.figure(figsize=(15, 10))
ax1 = fig.add_subplot(111)
ax1.set_xlabel("Time t")
ax1.set_title("relative energy error, rk4")

ax1.plot(times, relative_energy_error(t, all_y))
ax1.grid()

plt.savefig("p3_2.png")
