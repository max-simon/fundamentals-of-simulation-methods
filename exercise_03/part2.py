import numpy as np
import tree_alg as tree
from random import random
import sys


# number of particles
n = 100000

simulation = tree.Tree(1)

for i in range(0, n):
    success = simulation.insert_particle(np.array([-0.5 + random(), -0.5 + random(), -0.5 + random()]), 1/n)
    if not success:
        break


acc_old = simulation.calculate_acc_the_old_way(2)
acc_new = simulation.calculate_acc(2)

#simulation.print_tree()
#simulation.proof_all_particles()
print("old way: acc =", acc_old)
print("new way: acc =", acc_new)
print("eta = ", np.linalg.norm(acc_new - acc_old)/np.linalg.norm(acc_old))
