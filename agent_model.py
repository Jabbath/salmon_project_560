import numpy as np


class Salmon:
    def __init__(self, init_position, velocity, variance):
        """
        A salmon moving downstream in the river. A basic salmon starts
        as not infected by lice.
        :param init_position: The starting position of the salmon
        :param velocity: The average velocity with which the salmon moves downstream
        :param variance: The variance with which to add noise to the salmon's movements
        """
        self.position = init_position
        self.velocity = velocity
        self.variance = variance
        self.lice = []  # Attached lice

    def update_position(self, delta_t):
        """
        Move the salmon downstream with time difference delta_t. The
        salmon moves according to v*delta_t + B_{delta_t}, where B
        is a Brownian motion.
        :param delta_t: The time for which the salmon travels.
        :return: None
        """
        brownian_variance = (self.variance**2)*delta_t # See http://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-BM.pdf
        brownian_noise = np.random.normal(scale=brownian_variance)

        self.position = self.position + self.velocity*delta_t + brownian_noise


class River:
    def __init__(self, salmon_velocity, salmon_variance):
        '''
        Initializes the river. The river should contain salmon
        and lice. The initial distribution of lice should mirror
        the data given by Krkošek in his paper.
        :param salmon_velocity: The average velocity with which salmon
        move downstream.
        :param salmon_variance: The variance for a normally distributed
        noise term.
        '''
        self.salmon_velocity = salmon_velocity
        self.salmon_variance = salmon_variance
        self.salmon = []


def salmon_test():
    test_salmon = Salmon(0, 5, 1)
    delta_t = 1

    for i in range(10):
        test_salmon.update_position(delta_t)
        print(test_salmon.position)

salmon_test()

