import numpy as np
from numpy.linalg import norm


class Particle:
    def __init__(self, position, mass, tree):
        """
        Creates an particle.
        @param:
            position: numpy-array with x, y and z - position
            mass: mass of particle
        @return:
            Particle
        """
        self.position = position
        self.mass = mass

        self.global_particle_index = tree.register_particle(self)

        self.mother_node = 0

class Node:
    def __init__(self, center, length, tree, level, mother_node_index = -1):
        """
        Creates an particle.
        @param:
            center: numpy-array with x, y and z - position of center of node-block
            length: length of node-block
            tree: the Tree-Object the Node belongs to

        @return:
            Node
        """
        self.center = center
        self.length = length
        self.tree = tree

        self.com = None
        self.mass = 0
        self.subnodes = None
        self.mothernode_global_index = mother_node_index
        self.inner_particles = []

        self.level = level

        self.global_node_index = tree.register_node(self)
        # print("Registration of node %d" % self.global_node_index)


class Tree:
    def __init__(self, length, opening_threshold = 0.8):

        self.opening_threshold = opening_threshold

        self.all_nodes = []

        self.root = Node(np.array([0,0,0]), length, self, 0)

        self.all_particles = []

        self.multipoles_up_to_date = False

        self.max_nodes = 500000
        self.max_sublevels = 50
        self.softening = 1e-3

    # Register nodes and particles
    def register_node(self, node):
        self.all_nodes.append(node)
        return len(self.all_nodes) - 1
    def register_particle(self, particle):
        self.all_particles.append(particle)
        return len(self.all_particles) - 1

    def get_center_of_subnode(self, node, subnode_index):
        """
        Calculates the center of subnode with index subnode_index.
        @param:
            subnode_index: index of the subnode of interest
        @return:
            center: numpy-array with x, y and z of center
        """
        zfactor = 1 if subnode_index < 4 else -1
        yfactor = 1 if subnode_index%4 < 2 else -1
        xfactor = 1 if subnode_index%4 != 1 and subnode_index%4 != 2 else -1
        factors = np.array([xfactor, yfactor, zfactor])

        return node.center + factors*0.25*node.length

    def get_subnode_index(self, node, particle):
        """
        Return the index of the subnode the particle belongs to.
        @param:
            particle: particle to place
        @return:
            index: index of corresponding subnode
        """

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
        Creates subnodes
        @param
        @return
        """
        if len(self.all_nodes) + 8 > self.max_nodes: # if maximum of nodes is reached
            return False

        node.subnodes = []
        for i in range(8):
            new_node = Node(self.get_center_of_subnode(node, i), 0.5*node.length, self, node.level + 1, node.global_node_index)
            node.subnodes.append(new_node.global_node_index)

        return True

    def insert_particle(self, position, mass):
        # create a particle
        new_part = Particle(position, mass, self)

        # print("\nParticle %d created" % new_part.global_particle_index)

        # put particle into corresponding node
        current = self.root

        while len(current.inner_particles) > 0: # while current has at least one particle
            # print("Node %d has at least one particle (%d)" % (current.global_node_index, len(current.inner_particles)))

            if current.level >= self.max_sublevels: # if maximal sublevel is reached, put just on top of others in this node
                print("Maximal sublevel has reached. Particle %d was just dropped here without subnodes. Level = %d" % (new_part.global_particle_index, current.level))
                break

            if len(current.inner_particles) == 1: # if current has already one: create new subnodes and move existing particle to corresponding subnode
                success = self.create_subnodes(current)

                if not success:
                    print("Number of max nodes exceeded!")
                    self.all_particles.pop()
                    return False

                # print("Created new subnodes on node %d, subnode-level: %d (%d, %d, %d, %d, %d, %d, %d, %d)"%(current.global_node_index, current.level + 1, *current.subnodes))
                subnode_index = self.get_subnode_index(current, self.all_particles[current.inner_particles[0]])
                global_node_index = current.subnodes[subnode_index]
                self.all_nodes[global_node_index].inner_particles.append(current.inner_particles[0])
                # print("Moved particle %d to node %d" % (current.inner_particles[0], global_node_index))
                self.all_particles[current.inner_particles[0]].mother_node = global_node_index

            # from here there should be subnodes!
            current.inner_particles.append(new_part.global_particle_index)
            subnode_index = self.get_subnode_index(current, new_part)
            current = self.all_nodes[current.subnodes[subnode_index]]
            # print("Now look deeper into node %d" % current.global_node_index)

        current.inner_particles.append(new_part.global_particle_index)
        new_part.mother_node = current.global_node_index
        # print("Particle %d set into node %d" % (new_part.global_particle_index, current.global_node_index))

        self.multipoles_up_to_date = False
        return True


    def update_multipole_moments(self, node):
        """
        Updates the multipole moments of the node and subnodes recursively.
        @param
        @return
        """

        if node.subnodes != None: # calculate from subnodes
            # calculate multipole moments recursively
            for subnode_index in node.subnodes:
                subnode = self.all_nodes[subnode_index]
                self.update_multipole_moments(subnode)
            # calculate multipol moment from subnodes
            node.mass = 0
            node.com = np.array([0, 0, 0])
            # TODO: how to calculate the rest (quadrupoles)
            for subnode_index in node.subnodes:
                subnode = self.all_nodes[subnode_index]
                if subnode.mass == 0: # if subnode is empty: go to next subnode
                    continue
                node.mass += subnode.mass
                node.com = node.com + (subnode.mass * subnode.com)
            node.com = node.com/node.mass if node.mass != 0 else None

        else: # calculate from inner particles
            node.mass = 0
            node.com = np.array([0, 0, 0])
            # TODO: how to calculate the rest (quadrupole)
            for particle_index in node.inner_particles:
                particle = self.all_particles[particle_index]
                node.mass += particle.mass
                node.com = node.com + (particle.mass * particle.position)
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

    def get_acc_for(self, node, particle_index):

        # exclude empty nodes
        if node.mass == 0:
            return 0

        #load particle
        par1 = self.all_particles[particle_index]
        # get opening angle
        y = par1.position - node.com
        theta = self.get_opening_angle(node, y)

        if theta == None:
            return 0

        # or expression is, when maximum number of sublevels was reached
        if theta < self.opening_threshold or node.subnodes == None:
            acc = np.array([0, 0, 0])

            # monopole
            acc = acc - node.mass * y/((norm(y)**2 + self.softening**2)**(3/2))


            # quadrupole
            #q = self.calc_quadrupol_matrix(node)
            #acc = acc +  ((2 * np.dot(q, y) * np.linalg.norm(y)**5 + 10 * y * np.linalg.norm(y)**4 * np.dot(y, np.dot(q, y)))/(np.abs(y)**10))

        elif theta >= node.tree.opening_threshold:
            # open subnodes
            acc = np.array([0, 0, 0])

            for subnode_index in node.subnodes:
                acc = acc + self.get_acc_for(self.all_nodes[subnode_index], particle_index)

        else: # if something is going wrong :/
            acc = 0

        return acc


    def calculate_acc_the_old_way(self, particle_index):
        par1 = self.all_particles[particle_index]
        a = np.array([0, 0, 0])

        for particle in self.all_particles:
            if particle_index != particle.global_particle_index:
                a = a - particle.mass * (par1.position - particle.position)/((norm(particle.position - par1.position)**2 + self.softening**2)**(3/2))
        return a





    def calculate_acc(self, particle_index):
        if not self.multipoles_up_to_date:
            self.update_multipole_moments(self.root)
        return self.get_acc_for(self.root, particle_index)

    def proof_all_particles(self):
        for particle in self.all_particles:
            node = self.all_nodes[particle.mother_node]
            if (np.abs(node.center - particle.position) > np.array([node.length/2 for i in range(3)])).any():
                print("ERROR: Particle %d in wrong node" % particle.global_particle_index)
                print("   ", "Node-Center: (%1.6f, %1.6f, %1.6f), Node-Length: %1.6f, Particle-Position: (%1.6f, %1.6f, %1.6f)" % (*node.center, node.length, *particle.position))

    def print_node(self, node_index):
        # load node
        node = self.all_nodes[node_index]

        if node.mass == 0:
            return

        print("   "*node.level, " Node - %d p's, m = %1.8f (center: %1.6f, %1.6f, %1.6f, length = %1.6f)" % (len(node.inner_particles), node.mass, *node.center, node.length))

        if node.subnodes != None:
            print("   "*node.level, "   |")
            for subnode_index in node.subnodes:
                self.print_node(subnode_index)

        else:
            for particle_index in node.inner_particles:
                particle = self.all_particles[particle_index]
                print("   "*node.level, "   - Particle %d (%1.6f, %1.6f, %1.6f)" % (particle_index, *particle.position))


    def print_tree(self):
        self.print_node(0)
