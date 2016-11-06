import numpy as np
import tree_alg as tree
from random import random
import sys


# number of particles
n = 100

simulation = tree.Tree(1)

for i in range(0, n):
    success = simulation.insert_particle(np.array([-0.5 + random(), -0.5 + random(), -0.5 + random()]), 1/n)
    if not success:
        break

simulation.analyze(0.8)
