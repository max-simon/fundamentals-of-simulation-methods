import numpy as np
from matplotlib import pyplot as plt



def dtheta(t, theta):
    # theta is dimensionless temperature, tau dimensionless time
    # k includes all constants
    # see solution-sheet
    k = 0.2415
    if theta < 1:
        return -k*theta**10
    else:
        return -k*theta**(-0.5)

# converting values
def kelvin_to_theta(T):
    return T/20000
def theta_to_kelvin(theta):
    return theta*20000
def tau_to_time(tau):
    return tau*1e10
def time_to_tau(t):
    return t/1e10

# integrator
def rk2(f, theta, tau, dtau):
    k1 = f(tau, theta)
    k2 = f(tau + dtau, theta + dtau * k1)
    return theta + 0.5 * dtau * (k1 + k2)

# integrator with adaptive stepsize
# return used_dtau, theta, dtau_for_next_step
def adaptive_stepsize_integration(f, theta, tau, dtau, trunc_error=50, integrator = rk2, order_of_integrator = 2):
    truncation_error = kelvin_to_theta(trunc_error) # just for readability
    factor_for_very_much_less = 2**(order_of_integrator + 1)

    # large step
    theta_1 = integrator(f, theta, tau, dtau)
    # two smaller step
    theta_tmp = integrator(f, theta, tau, dtau/2)
    theta_2 = integrator(f, theta_tmp, tau, dtau/2)

    if np.abs(theta_2 - theta_1) < truncation_error/factor_for_very_much_less: # if difference is smaller than error
        return dtau, theta_2, dtau*2 # adjust dtau and give more precise theta_2 back
    elif np.abs(theta_2 - theta_1) <= truncation_error:
        return dtau, theta_2, dtau # accept step, but do not change dtau
    else: # if error is to big
        return adaptive_stepsize_integration(f, theta, tau, dtau/2, trunc_error, integrator, order_of_integrator) # try again with smaller stepsize



###################
#### 2a and b #####
###################

print("----- Exercise 2a and 2b -----")

# create list of dtau's to check
dtaus = [time_to_tau(i*1e10) for i in range(1, 10)]

# initialise plot figure
fig, (ax1) = plt.subplots(1, 1, figsize=(20, 10))
ax1 = fig.add_subplot(111)

for_exercise_part_c_tau = np.array(0)
for_exercise_part_c_theta = np.array(kelvin_to_theta(1e7))
for_exercise_part_c_is = 0

for dtau in dtaus:
    # Presets
    tau = 0
    theta = kelvin_to_theta(1e7)

    all_theta = np.array([theta])
    all_tau = np.array([tau])
    integration_steps = 0

    while theta > 0.3: # 6000K in theta: 0.3
        tau += dtau
        theta = rk2(dtheta, theta, tau, dtau)

        # this loop is just because I do not want to integrate again in part c.
        if dtau == 1:
            for_exercise_part_c_tau = np.append(for_exercise_part_c_tau, tau)
            for_exercise_part_c_theta = np.append(for_exercise_part_c_theta, theta)

        all_theta = np.append(all_theta, theta)
        all_tau = np.append(all_tau, tau)

        integration_steps += 1

    if dtau == 1:
        for_exercise_part_c_is = integration_steps

    print("dtau =", dtau, " - Integration finished.\n\tIntegration steps:", integration_steps)
    ax1.plot(all_tau, all_theta, label="Simulation with RK2\n$\\Delta\\tau$ = {:1.0f}\nIteration steps: {:d}".format(dtau, integration_steps))
    print("\tEvolution plotted.")


ax1.set_yscale("log")
ax1.grid()
ax1.set_xlabel("Time in $10^{10}$sec")
ax1.set_ylabel("Temperatur in T$_0$")
ax1.set_title("Time evolution of temperature")
ax1.legend()
plt.savefig("p2_1.png")



###################
####    2c    #####
###################

print("----- Exercise 2c -----")

# Presets
tau = 0
dtau_next_step = time_to_tau(4e10)
theta = kelvin_to_theta(1e7)

all_theta = np.array([theta])
all_tau = np.array([tau])
all_dtau = np.array([dtau_next_step])
integration_steps = 0

while theta > 0.3: # 6000K in theta: 0.3
    dtau, theta, dtau_next_step = adaptive_stepsize_integration(dtheta, theta, tau, dtau_next_step, 50, rk2, 2)
    tau += dtau

    all_theta = np.append(all_theta, theta)
    all_tau = np.append(all_tau, tau)
    all_dtau = np.append(all_dtau, dtau)

    integration_steps += 1

print("Integration with adaptive stepsize finished.\n\tIntegration steps:", integration_steps)


# initialise plot figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), sharex=True)

ax1.plot(all_tau, all_theta, label="Integration with adaptive stepsize\nIteration steps: {:d}".format(integration_steps))
ax1.plot(for_exercise_part_c_tau, for_exercise_part_c_theta, label="Simulation with RK2\n$\\Delta\\tau$ = {:1.0f}\nIteration steps: {:d}".format(1, for_exercise_part_c_is))
ax1.set_yscale("log")
ax1.grid()
ax1.set_xlabel("Time in $10^{10}$sec")
ax1.set_ylabel("Temperatur in T$_0$")
ax1.set_title("Time evolution of temperature")
ax1.legend()

# show stepsize evolution
ax2.plot(all_tau, all_dtau)
ax2.set_yscale("log")
ax2.set_xlabel("Time in $10^{10}$sec")
ax2.set_ylabel("Stepsize $\\Delta \\tau$")

plt.savefig("p2_2.png")

print("\tEvolution plotted.")
