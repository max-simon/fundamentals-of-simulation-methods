import numpy as np
from numpy.linalg import norm
from time import time

import sys


class Particle:
    def __init__(self, position, mass, tree):
        """
        Creates a particle.
        @param:
            position: numpy-array with x, y and z - position
            mass: mass of particle
            tree: the Tree-Object the particle belongs to
        @return:
            Particle
        """
        self.position = position # np-array with position of particle
        self.mass = mass # mass of particle

        self.global_particle_index = tree.register_particle(self) # global index referring to this particle

        self.mother_node = 0 # index of closest node the particle belongs to

class Node:
    def __init__(self, center, length, tree, level, mother_node_index = -1):
        """
        Creates a node.
        @param:
            center: numpy-array with x, y and z - position of center of node-block
            length: length of node-block
            tree: the Tree-Object the Node belongs to
            level: sublevel of this node
            mother_node_index: index of nearest node

        @return:
            Node
        """
        self.center = center # np-array of center of the Node
        self.length = length # length of node
        self.tree = tree # tree, to which the node belongs

        self.com = None # np-array of center of mass, will be calculated
        self.mass = 0 # mass, will be calculated
        self.quadrupole = 0

        self.subnodes = None # later a list of indices of the subnodes
        self.mothernode_global_index = mother_node_index # index of nearest node
        self.inner_particles = [] # list of particles, which are inside this node (all!)

        self.level = level # sublevel of this node

        self.global_node_index = tree.register_node(self) # global index referring to this node
        # print("Registration of node %d" % self.global_node_index)


class Tree:
    def __init__(self, init_length):
        """
        Creates a tree.
        @param:
            init_length: length of root-node

        @return:
            Tree
        """
        self.opening_threshold = 0
        self.all_nodes = [] # contains all nodes

        self.root = Node(np.array([0,0,0]), init_length, self, 0) # create root-node

        self.all_particles = [] # contains all particles
        self.multipoles_up_to_date = False # a boolean, are multipoles up to date?

        self.max_nodes = 500000 # maximal nodes
        self.max_sublevels = 50 # maximal depth
        self.softening = 1e-3 # softening factor

    def register_node(self, node):
        """
        Register a node in the tree.
        @param:
            node: node to register

        @return:
            unique index of node
        """
        self.all_nodes.append(node)
        return len(self.all_nodes) - 1

    def register_particle(self, particle):
        """
        Register a particle in the tree.
        @param:
            particle: particle to register

        @return:
            unique index of particle
        """
        self.all_particles.append(particle)
        return len(self.all_particles) - 1

    def get_center_of_subnode(self, node, subnode_index):
        """
        Calculates the center of subnode with index subnode_index.
        @param:
            node: node for which the subnode is
            subnode_index: index of the subnode of interest
        @return:
            center: numpy-array with x, y and z of center
        """
        # create an array of signs
        zfactor = 1 if subnode_index < 4 else -1
        yfactor = 1 if subnode_index%4 < 2 else -1
        xfactor = 1 if subnode_index%4 != 1 and subnode_index%4 != 2 else -1
        factors = np.array([xfactor, yfactor, zfactor])

        return node.center + factors*0.25*node.length

    def get_subnode_index(self, node, particle):
        """
        Get the index of subnode the particle should be in.
        @param:
            node: node the subnodes belong to
            particle: particle to place
        @return:
            index: return the index of the corresponding subnode
        """
        # Trust me :D
        index = 0
        if particle.position[0] > node.center[0] and particle.position[1] < node.center[1]:
            index = 3
        if particle.position[0] < node.center[0]:
            index = 1
            if particle.position[1] < node.center[1]:
                index = 2
        if particle.position[2] < node.center[2]:
            index += 4
        return index

    def create_subnodes(self, node):
        """
        Creates subnodes in node, if it is possible according to maximal nodes number.
        @param
            node: the node the subnodes should belong to
        @return
            bool: subnodes have been added or not
        """
        if len(self.all_nodes) + 8 > self.max_nodes: # if maximum of nodes is reached
            return False # do not make subnodes

        node.subnodes = [] # initialize subnodes
        for i in range(8): # create 8 new subnodes and append them to the node
            new_node = Node(self.get_center_of_subnode(node, i), 0.5*node.length, self, node.level + 1, node.global_node_index)
            node.subnodes.append(new_node.global_node_index)
        return True

    def insert_particle(self, position, mass):
        """
        Insert a particle into the tree.
        @param
            position: position of the particle
            mass: mass of particle
        @return
            bool: particle has been added or not
        """
        # create a particle
        new_part = Particle(position, mass, self)

        # print("\nParticle %d created" % new_part.global_particle_index)

        current = self.root # initialize following loop
        while len(current.inner_particles) > 0: # while current has at least one particle
            # print("Node %d has at least one particle (%d)" % (current.global_node_index, len(current.inner_particles)))

            if current.level >= self.max_sublevels: # if maximal sublevel is reached, put just on top of others in this node
                print("Maximal sublevel has reached. Particle %d was just dropped here without subnodes. Level = %d" % (new_part.global_particle_index, current.level))
                break # by just breaking the loop

            if len(current.inner_particles) == 1: # if current has already one: create new subnodes and move existing particle to corresponding subnode
                success = self.create_subnodes(current)

                if not success:
                    print("Number of max nodes exceeded!")
                    self.all_particles.pop() # remove particle
                    return False

                # print("Created new subnodes on node %d, subnode-level: %d (%d, %d, %d, %d, %d, %d, %d, %d)"%(current.global_node_index, current.level + 1, *current.subnodes))

                # move existing particle to the corresponding subnode
                subnode_index = self.get_subnode_index(current, self.all_particles[current.inner_particles[0]])
                global_node_index = current.subnodes[subnode_index] # get index
                self.all_nodes[global_node_index].inner_particles.append(current.inner_particles[0]) # move particle
                self.all_particles[current.inner_particles[0]].mother_node = global_node_index # say particle, that it has been moved
                # print("Moved particle %d to node %d" % (current.inner_particles[0], global_node_index))

            # from here there should be subnodes!

            # move new particle to corresponding subnode
            current.inner_particles.append(new_part.global_particle_index)

            # set new current
            subnode_index = self.get_subnode_index(current, new_part)
            current = self.all_nodes[current.subnodes[subnode_index]]
            # print("Now look deeper into node %d" % current.global_node_index)

        # after the loop: we reached a blank node

        # so put new particle to subnode
        current.inner_particles.append(new_part.global_particle_index) # move particle
        new_part.mother_node = current.global_node_index # say particle that it has been moved
        # print("Particle %d set into node %d" % (new_part.global_particle_index, current.global_node_index))

        self.multipoles_up_to_date = False
        return True

    def update_multipole_moments(self, node):
        """
        Updates the multipole moments of the node and subnodes recursively.
        @param:
            node: the starting node
        @return
        """
        if node.subnodes != None: # if subnodes exists: create from subnodes

            # calculate multipole moments recursively
            for subnode_index in node.subnodes:
                subnode = self.all_nodes[subnode_index]
                self.update_multipole_moments(subnode)

            # calculate multipol moment from subnodes
            node.mass = 0
            node.com = np.array([0, 0, 0])

            # mass and com
            for subnode_index in node.subnodes:
                subnode = self.all_nodes[subnode_index]
                if subnode.mass == 0: # if subnode is empty: go to next subnode
                    continue
                node.mass += subnode.mass
                node.com = node.com + (subnode.mass * subnode.com)

            # TODO: how to calculate the rest (quadrupoles)
            # now uses all inner particles.
            node.quadrupole = self.calc_quadrupol_matrix(node)

            # only com if mass is not 0
            node.com = node.com/node.mass if node.mass != 0 else None

        else: # calculate from inner particles
            node.mass = 0
            node.com = np.array([0, 0, 0])

            # mass and com
            for particle_index in node.inner_particles:
                particle = self.all_particles[particle_index]
                node.mass += particle.mass
                node.com = node.com + (particle.mass * particle.position)

            # TODO: how to calculate the rest (quadrupoles)
            # now uses all inner particles.
            node.quadrupole = self.calc_quadrupol_matrix(node)

            # only com if mass is not 0
            node.com = node.com/node.mass if node.mass != 0 else None

    def get_opening_angle(self, node, y):
        """
        Calculates the opening angle according to pos. First the update_com_and_mass has to be called!
        @param
            pos: numpy-array with x, y and z giving the position the angle corresponds to
        @return
            angle: opening-angle in rad
        """
        return node.length / norm(y) if norm(y) != 0 else None

    # maybe useless
    def calc_quadrupol_matrix(self, node):
        """
        Calculates the quadrupole-matrix.
        @param
        @return
            matrix: numpy-matrix
        """
        q_m = np.array([])
        for i in range(3):
            for j in range(i, 3):
                q = 0
                # loop over all inner particles -> only over subnodes?
                for particle_index in node.inner_particles:
                    particle = self.all_particles[particle_index]
                    cron_delta = 1 if i == j else 0
                    q += particle.mass*(3*(np.dot(node.com[i] - particle.position[i], node.com[j] - particle.position[j]) - cron_delta*norm(node.com - particle.position)**2))
                q_m = np.append(q_m, q)

        return np.array([[q_m[0], q_m[1], q_m[2]], [q_m[1], q_m[3], q_m[4]], [q_m[2], q_m[4], q_m[5]]])

    def get_acc_for(self, node, particle_index, term_counter = 0):
        """
        Calculates the acceleration on a particle starting at node (recursively).
        @param
            node: node to start from
            particle_index: particle on which the acceleration appear
            term_counter: counts how many nodes were used
        @return
            acc: the acceleration
        """
        # exclude empty nodes
        if node.mass == 0:
            return 0, term_counter

        #load particle
        par1 = self.all_particles[particle_index]

        # get opening angle
        y = par1.position - node.com
        theta = self.get_opening_angle(node, y)

        # if y = 0, excluded because infinity
        if theta == None:
            return 0, term_counter

        # or expression is, when maximum number of sublevels was reached
        # and there are several particles in the node
        if theta < self.opening_threshold or node.subnodes == None:

            acc = np.array([0, 0, 0])

            # monopole
            acc = acc - node.mass * y/((norm(y)**2 + self.softening**2)**(3/2))

            # quadrupole
            #acc = acc +  ((2 * np.dot(node.quadrupole, y) * norm(y)**5 + 10 * y * norm(y)**4 * np.dot(y, np.dot(node.quadrupole, y)))/(norm(y)**10))

            term_counter += 1

        elif theta >= node.tree.opening_threshold: # open subnodes
            acc = np.array([0, 0, 0])

            for subnode_index in node.subnodes:
                new_acc, term_counter = self.get_acc_for(self.all_nodes[subnode_index], particle_index, term_counter)
                acc = acc + new_acc

        else: # this should never happen :D
            acc = 0

        return acc, term_counter

    def calculate_acc_exact(self, particle_index):
        """
        Calculate the exact acceleration by summing up all.
        @param:
            particle_index: the index of the particle the acceleration depends on.
        @return:
            acc: the acceleration
        """

        par1 = self.all_particles[particle_index] # load the particle

        acc = np.array([0, 0, 0])

        for particle in self.all_particles:
            if particle_index != particle.global_particle_index: # if not is the same particle
                # increase acceleration
                acc = acc - particle.mass * (par1.position - particle.position)/((norm(particle.position - par1.position)**2 + self.softening**2)**(3/2))
        return acc

    def calculate_acc(self, particle_index, opening_threshold):
        """
        Calculate the acceleration by tree method.
        @param:
            particle_index: the index of the particle the acceleration depends on.
            opening_threshold: the threshold used for decide if subnodes should be opened or not
        @return:
            acc: the acceleration
        """
        self.opening_threshold = opening_threshold
        if not self.multipoles_up_to_date: # are multipoles already calculated
            self.update_multipole_moments(self.root)
            self.multipoles_up_to_date = True
        # calculate by calling get_acc_for
        acc, terms_used = self.get_acc_for(self.root, particle_index, 0)
        return acc, terms_used

    def proof_all_particles(self):
        """
        Proofs, if all particles are in the right subnode.
        @param:
        @return:
        """
        for particle in self.all_particles:
            # load node the particle claims to belong to
            node = self.all_nodes[particle.mother_node]
            # if the particle is out of this node
            if (np.abs(node.center - particle.position) > np.array([node.length/2 for i in range(3)])).any():
                print("ERROR: Particle %d in wrong node" % particle.global_particle_index)
                print("   ", "Node-Center: (%1.6f, %1.6f, %1.6f), Node-Length: %1.6f, Particle-Position: (%1.6f, %1.6f, %1.6f)" % (*node.center, node.length, *particle.position))
            # if there are subnodes the particle can also belong to
            if node.subnodes != None:
                print("WARNING: The node %d has subnodes. The particle %d may should be in one of those." % (particle.mother_node, particle.global_particle_index))

    def print_node(self, node_index):
        """
        Prints the tree from given node (recursively)
        @param:
            node_index: the node from which the tree should be printed
        @return:
        """
        # load node
        node = self.all_nodes[node_index]
        # skip all empty nodes
        if node.mass == 0:
            return

        print("   "*node.level, " Node - %d p's, m = %1.8f (center: %1.6f, %1.6f, %1.6f, length = %1.6f)" % (len(node.inner_particles), node.mass, *node.center, node.length))

        # if node have subnodes: print them
        if node.subnodes != None:
            print("   "*node.level, "   |")
            for subnode_index in node.subnodes:
                self.print_node(subnode_index)

        else:
            # print all inner particles
            for particle_index in node.inner_particles:
                particle = self.all_particles[particle_index]
                print("   "*node.level, "   - Particle %d (%1.6f, %1.6f, %1.6f)" % (particle_index, *particle.position))

    def print_tree(self):
        """
        Prints the complete tree
        @param:
        @return:
        """
        self.print_node(0)

    def analyze(self, threshold, single_particles = False, proof_particles = False):
        """
        Analysing the tree algorithm in comparism to exact method.
        @param:
            threshold: the threshold for the opening angle
            single_particles: list of indices particles one want to analyze. If False, all particles are used for analysation
        """
        # initialise list of particles used for the analysation
        if single_particles != False:
            particles = single_particles
        else:
            particles = range(len(self.all_particles))

        # proof particles?
        if proof_particles:
            print("Proof positions of all particles...")
            self.proof_all_particles()

        # exact
        print("Started exact calculation...")
        all_exact = []
        t0 = time() # take starting time point
        counter = 0
        for particle_index in particles: # loop over particles in the list
            progressBar("\tProgress", counter, len(particles)) # update command line
            acc_exact = self.calculate_acc_exact(particle_index)
            all_exact.append(acc_exact) # save for later analysation
            counter += 1
        t1 = time() # take end time point
        time_of_exact = t1 - t0
        print("\tFinished.")

        # tree
        print("Started calculation with tree method...")
        total_terms = 0
        all_tree = []
        t0 = time()
        counter = 0
        for particle_index in particles: # loop over particles in the list
            progressBar("\tProgress", counter, len(particles)) # update command line
            acc_tree, terms_used = self.calculate_acc(particle_index, threshold)
            all_tree.append(acc_tree) # save for later analysation
            total_terms += terms_used # save for later analysation
            counter += 1
        t1 = time()
        time_of_tree = t1 - t0
        total_terms = total_terms/len(particles)
        print("\tFinished")

        # Analysation
        eta = 0
        for i in range(len(particles)):
            eta += norm(all_exact[i] - all_tree[i])/norm(all_exact[i]) # calculate mean error
        eta = eta/len(particles)

        # output
        print("Analysation: N = {:d} particles, threshold = {:1.3f}, total mass: {:3.3f}".format(len(particles), self.opening_threshold, self.root.mass))
        print("\tExact:")
        print("\t\ttime: {:1.6f}".format(time_of_exact))
        print("\tTree:")
        print("\t\ttime: {:1.6f}".format(time_of_tree))
        print("\t\tmean relative error: {:2.5f}".format(eta))
        print("\t\tmean used nodes: {:4.2f}".format(total_terms))

        return len(particles), threshold, time_of_exact, time_of_tree, eta 

    def reset(self):
        """
        Resets the complete tree.
        @param
        @return
        """
        # reset to initial values
        self.opening_threshold = 0
        self.all_nodes = []
        self.root = None
        self.all_particles = []
        self.multipoles_up_to_date = False
        self.max_nodes = 500000
        self.max_sublevels = 50
        self.softening = 1e-3

    @staticmethod
    def init_a_tree(n, distribution, mass_function, init_length = 1, tree = None):
        """
        Initialise a tree with n particles and a given distribution function.
        @param:
            n: number of particles
            distribution: probability-distribution for placement of the particles
            mass of particles: function of n, x, y and z, which return the mass of the particle
        @return
        """
        # create a tree object
        if tree == None:
            tree = Tree(init_length)
            print("Created a new tree...")

        t0 = time()
        for i in range(n):
            pos = np.array([-init_length/2 + distribution()*init_length, -init_length/2 + distribution()*init_length, -init_length/2 + distribution()*init_length])
            success = tree.insert_particle(pos, mass_function(n, *pos))
            if not success:
                break
        t1 = time()

        print("\tInserted {:d} particles".format(len(tree.all_particles)))
        print("\t{:d} nodes created".format(len(tree.all_nodes)))
        print("\tTime needed: {:3.3f}".format(t1 - t0))

        return tree


def progressBar(title, value, endvalue, bar_length=20):
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\r{0}: [{1}] {2}%".format(title, arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()
