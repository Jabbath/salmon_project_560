import numpy as np
import scipy

class Salmon:
    def __init__(self, init_position, velocity, sigma):
        """
        A salmon moving downstream in the river. A basic salmon starts
        as not infected by lice.
        :param init_position: The starting position of the salmon
        :param velocity: The average velocity with which the salmon moves downstream
        :param variance: The variance with which to add noise to the salmon's movements
        """
        self.position = init_position
        self.velocity = velocity
        self.variance = sigma
        self.lice = []  # Attached lice

    def update_position(self, delta_t):
        """
        Move the salmon downstream with time difference delta_t. The
        salmon moves according to v*delta_t + B_{delta_t}, where B
        is a Brownian motion.
        :param delta_t: The time for which the salmon travels.
        :return: None
        """
        brownian_variance = (self.sigma**2)*delta_t # See http://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-BM.pdf
        brownian_noise = np.random.normal(scale=brownian_variance)

        self.position = self.position + self.velocity*delta_t + brownian_noise


class River:
    def __init__(self, salmon_velocity, salmon_sigma, louse_velocity, louse_sigma,
                 salmon_spawn_rate, louse_farm_rate, louse_ambient_rate, infection_lambda,
                 start_x=-20, end_x=20):
        '''
        Initializes the river. The river should contain salmon
        and lice. The initial distribution of lice should mirror
        the data given by Krkošek in his paper.
        :param salmon_velocity: The average velocity with which salmon
        move downstream.
        :param salmon_variance: The sigma for a normally distributed
        noise term.
        '''
        self.salmon = []
        self.salmon_velocity = salmon_velocity
        self.salmon_sigma = salmon_sigma

        self.lice = []
        self.louse_velocity = louse_velocity
        self.louse_sigma = louse_sigma

        self.salmon_spawn_rate = salmon_spawn_rate
        self.louse_farm_rate = louse_farm_rate
        self.louse_ambient_rate = louse_ambient_rate

        self.infection_lambda = infection_lambda

        self.start_x = start_x
        self.end_x = end_x

    def update(self, delta_t):
        """
        Updates all elements of the river. First, we spawn new salmon and lice.
        Then we update the positions of the salmon, and lice. We then age the lice.
        Finally, we compute new infections.
        :param delta_t: The timestep for which we update
        :return:
        """
        # Spawn salmon
        num_salmon_to_spawn = int(self.salmon_spawn_rate*delta_t)

        for i in range(num_salmon_to_spawn):
            new_salmon = Salmon(self.start_x, self.salmon_velocity, self.salmon_sigma)
            self.salmon.append(new_salmon)

        # Spawn farm lice
        # num_farm_lice_to_spawn = int(self.louse_farm_rate*delta_t)
        #
        # for i in range(num_farm_lice_to_spawn):
        #     new_louse = Louse(0, self.louse_velocity, self.louse_sigma)
        #     self.lice.append(new_louse)

        # Update salmon positions
        for salmon in self.salmon:
            salmon.update_position(delta_t)

        # Update lice positions
        # for louse in self.lice:
        #     louse.update_position(delta_t)

        # Age lice
        # for louse in self.lice:
        #     louse.age(delta_t)

        # Calculate infections
        #self.update_infections(delta_t)

    def update_infections(self, delta_t):
        """
        Computes new infections between lice and salmon. For each salmon within a
        fixed distance of a louse, we draw from an exp. dist. for a time. If the
        time is less than the timestep size, we infect.
        :param delta_t:
        :return:
        """
        # Calculate a louse-salmon distance matrix
        louse_positions = np.array([louse.position for louse in self.lice], shape=(len(self.lice), 1))
        salmon_positions = np.array([salmon.position for salmon in self.salmon], shape=(len(self.salmon), 1))

        dist_matrix = scipy.spatial.distance_matrix(louse_positions, salmon_positions)
        dist_matrix = np.abs(dist_matrix)  # We only care about the magnitude of distances

        # Find which salmon are within the distance threshold
        threshold = 3*self.salmon_sigma*np.sqrt(delta_t)
        threshold_matrix = dist_matrix < threshold

        # For each louse and each salmon within the threshold, check if we infect
        for louse_idx in range(threshold_matrix.shape[0]):
            possible_targets = np.where(threshold_matrix[louse_idx, :])[0]

            # Sample from the exp. dist to check if we infect
            for target_salmon_idx in possible_targets:
                sampled_time = np.random.exponential(scale=1/self.infection_lambda)

                # If the salmon is infected, attach the louse to it
                if sampled_time <= delta_t:
                    salmon = self.salmon[target_salmon_idx]
                    louse = self.lice[louse_idx]

                    louse.attach(salmon)
                    break

class Lice:
    def __init__(self, init_position, velocity_lice, sigma_lice, parent, attached, age, stage):
        """
        A salmon moving downstream in the river. A basic salmon starts
        as not infected by lice.
        :param init_position: The starting position of the lice
        :param velocity: The average velocity with which the lice moves downstream
        :param variance: The variance with which to add noise to the lice's movements
        :param parent: The infected salmon
        :param attached: True or False
        :param age: age of louse
        :param stage: stage of life of louse
        """
        self.position = init_position
        self.velocity = velocity_lice
        self.variance = sigma_lice
        self.parent = parent
        self.attached = attached
        self.age = age
        self.stage = stage

    def update_position(self, delta_t):
        """
        Move the lice downstream with time difference delta_t. The
        lice moves according to v*delta_t + B_{delta_t}, where B
        is a Brownian motion.
        :param delta_t: The time for which the salmon travels.
        :return: None
        """
        brownian_variance = (self.sigma**2)*delta_t # See http://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-BM.pdf
        brownian_noise = np.random.normal(scale=brownian_variance)
        if 
        self.position = self.position + self.velocity*delta_t + brownian_noise
    def add_parent(Salmon):
        """
        Add a parent salmon
        """
       self.parent = Salmon
        


def salmon_test():
    test_salmon = Salmon(0, 5, 2)
    delta_t = 1

    for i in range(10):
        test_salmon.update_position(delta_t)
        print(test_salmon.position)

salmon_test()


