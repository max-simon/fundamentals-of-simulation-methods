import numpy as np
import tree_alg as tree
from random import random

from matplotlib import pyplot as plt


# number of particles
ns = [200, 300, 400] # n's we want to test
thetas = [0.4] # thetas we want to test

probability_distribution = random # use a uniform distribution for placements of the particles
mass_function = lambda n, x, y, z: 1/n # every particle should have the same mass of 1/n

# variables for analysation
all_n = np.array([])
all_time_exact = np.array([])
all_time_tree = np.array([])
all_eta = np.array([])
all_theta = np.array([])

# loop over grid
for n in ns:
    for theta in thetas:
        print("\n----- N = {:d}, Theta = {:1.2f} -----".format(n, theta))
        simulation = tree.Tree.init_a_tree(n, probability_distribution, mass_function)
        n, threshold, time_of_exact, time_of_tree, eta = simulation.analyze(theta)

        all_n = np.append(all_n, n)
        all_theta = np.append(all_theta, theta)
        all_time_exact = np.append(all_time_exact, time_of_exact)
        all_time_tree = np.append(all_time_tree, time_of_tree)
        all_eta = np.append(all_eta, eta)


#### Plotting

theta_of_interest = 0.4
print("\n\nPlot evaluation for Theta = {:1.2f}".format(theta_of_interest))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), sharex=True)

ax1.plot(all_n[all_theta == theta_of_interest], all_time_exact[all_theta == theta_of_interest], label="Time of exact calculation")
ax1.plot(all_n[all_theta == theta_of_interest], all_time_tree[all_theta == theta_of_interest], label="Time of tree calculation")
ax1.set_ylabel("seconds")
ax1.set_title("Comparism of durations")
ax1.legend(loc=2)
ax2.plot(all_n[all_theta == theta_of_interest], all_eta[all_theta == theta_of_interest], label="Mean relative error")
ax2.legend(loc=2)
ax2.set_title("Relative mean error")

plt.savefig("evaluation.png")


#TODO: Welche Skala?

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), sharex=True)

ax1.plot(all_n[all_theta == theta_of_interest], all_time_exact[all_theta == theta_of_interest], label="Time of exact calculation")
ax1.plot(all_n[all_theta == theta_of_interest], all_time_tree[all_theta == theta_of_interest], label="Time of tree calculation")
ax1.set_ylabel("seconds")
ax1.set_title("Comparism of durations")
ax1.set_xscale("log")
ax1.legend(loc=2)
ax2.plot(all_n[all_theta == theta_of_interest], all_eta[all_theta == theta_of_interest], label="Mean relative error")
ax2.legend(loc=2)
ax2.set_title("Relative mean error")

plt.savefig("evaluation1.png")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), sharex=True)

ax1.plot(all_n[all_theta == theta_of_interest], all_time_exact[all_theta == theta_of_interest], label="Time of exact calculation")
ax1.plot(all_n[all_theta == theta_of_interest], all_time_tree[all_theta == theta_of_interest], label="Time of tree calculation")
ax1.set_ylabel("seconds")
ax1.set_title("Comparism of durations")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.legend(loc=2)
ax2.plot(all_n[all_theta == theta_of_interest], all_eta[all_theta == theta_of_interest], label="Mean relative error")
ax2.legend(loc=2)
ax2.set_title("Relative mean error")

plt.savefig("evaluation2.png")
