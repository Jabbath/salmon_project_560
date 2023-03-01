import numpy as np
import scipy
from scipy import spatial
from tqdm import tqdm
import matplotlib.pyplot as plt

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
        self.infected = False

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

class Louse:
    def __init__(self, init_position, velocity_lice, variance_lice, parent=None, attached=False, age=0, stage="copepodid",
                 alive=True, gave_birth= False, offspring_number=20):
        """
        A salmon moving downstream in the river. A basic salmon starts
        as not infected by lice.
        :param init_position: The starting position of the lice
        :param velocity: The average velocity with which the lice moves downstream
        :param variance: The variance with which to add noise to the lice's movements
        :param parent: The infected salmon
        :param attached: True or False
        :param age: age of louse
        :param stage: stage of life of louse: copepodid, chalimus, adulting
        """
        self.position = init_position
        self.velocity = velocity_lice
        self.variance = variance_lice
        self.parent = parent
        self.attached = attached
        self.age = age
        self.stage = stage
        self.alive = alive
        self.gave_birth = gave_birth
        self.offspring_number = offspring_number

    def update_position(self, delta_t):
        """
        Move the lice downstream with time difference delta_t. The
        lice moves according to v*delta_t + B_{delta_t}, where B
        is a Brownian motion when it is not attached. Otherwise,
        its position is the position of the salmon it is attached to.

        :param delta_t: The time for which the salmon travels.
        :return: None
        """
        brownian_variance = (self.variance**2)*delta_t # See http://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-BM.pdf
        brownian_noise = np.random.normal(scale=brownian_variance)
        if not self.attached:
            self.position = self.position + self.velocity*delta_t + brownian_noise
        elif self.attached:
            self.position = self.parent.position
    def attach(self, salmon):
        """
        Add a parent salmon.

        :param salmon: The salmon to which the louse attaches
        """
        self.parent = salmon
        salmon.lice.append(self)
        salmon.infected = True
        self.attached = True
        self.age = 0

    def increase_age(self, delta_t):
        """
        Ages the salmon and moves it between life stages.

        :param delta_t: The time step to age by (in days)
        :return:
        """

        # If the salmon doesn't become attached in a set amount of days, it dies
        if not self.attached:
            if self.age > 10:
                self.alive = False
            else :
                self.age += delta_t

        # Otherwise age the louse
        elif self.attached:
            self.age += delta_t

            # Age the louse and kill it with probability according the Krkosek supplement
            if (10 < self.age and self.age <= 35) and (self.stage == 'copepodid'):
                # Draw from uniform to see if we survive
                survival_prob = 0.8804
                draw = np.random.uniform()

                if draw > survival_prob:
                    self.alive = False

                self.stage = 'chalimus'
            elif (35 < self.age) and (self.stage == 'chalimus'):
                # Draw from uniform to see if we survive
                survival_prob = 0.2586
                draw = np.random.uniform()

                if draw > survival_prob:
                    self.alive = False

                self.stage = 'adulting'

    def give_birth(self):
        """
        Allows the louse to create new lice at it's current location. Can only be
        used once.

        :return: new_lice_array: An array of the new lice that this louse gave birth to
        """
        if not self.gave_birth:
            new_lice_array = []

            for i in range(self.offspring_number):
                new_louse = Louse(self.position, self.velocity, self.variance)
                new_lice_array.append(new_louse)

            self.gave_birth = True

            return new_lice_array
        else:
            return []

class River:
    def __init__(self, salmon_velocity, salmon_sigma, louse_velocity, louse_sigma,
                 salmon_spawn_rate, louse_farm_rate, louse_ambient_rate, infection_lambda,
                 offspring_number=20, start_x=-20, end_x=20):
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
        self.louse_offspring = offspring_number

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
        expected_num_salmon = self.salmon_spawn_rate*delta_t
        num_salmon_to_spawn = np.random.poisson(lam=expected_num_salmon)

        for i in range(num_salmon_to_spawn):
            new_salmon = Salmon(self.start_x, self.salmon_velocity, self.salmon_sigma)
            self.salmon.append(new_salmon)

        # Spawn farm lice
        expected_num_farm_lice = self.louse_farm_rate*delta_t
        num_farm_lice_to_spawn = np.random.poisson(lam=expected_num_farm_lice)

        for i in range(num_farm_lice_to_spawn):
            new_louse = Louse(0, self.louse_velocity, self.louse_sigma)
            self.lice.append(new_louse)

        # Spawn wild lice
        expected_num_wild_lice = self.louse_ambient_rate*delta_t
        num_wild_lice_to_spawn = np.random.poisson(lam=expected_num_wild_lice)

        for i in range(num_wild_lice_to_spawn):
            location = np.random.uniform(self.start_x, self.end_x)
            new_louse = Louse(location, self.louse_velocity, self.louse_sigma)
            self.lice.append(new_louse)

        # Update salmon positions
        for salmon in self.salmon:
            salmon.update_position(delta_t)

        #Update lice positions
        for louse in self.lice:
            louse.update_position(delta_t)

        # Age lice
        for louse in self.lice:
            louse.increase_age(delta_t)

            # Check if the louse is old enough to give birth and make new lice if it is
            if louse.age >= 57 and not louse.gave_birth:
                new_lice = louse.give_birth()
                self.lice.extend(new_lice)

        # Subset the lice to only those that are still alive
        self.lice = [louse for louse in self.lice if louse.alive]

        # Remove any dead lice that might be attached to salmon
        for salmon in self.salmon:
            if salmon.infected:
                new_lice = [louse for louse in salmon.lice if louse.alive]
                salmon.lice = new_lice

                if len(new_lice) == 0:
                    salmon.infected = False

        # Subset salmon and lice to only those within the viewing area
        self.salmon = [salmon for salmon in self.salmon if salmon.position <= self.end_x]
        self.lice = [louse for louse in self.lice if louse.position <= self.end_x]

        # Calculate infections
        self.update_infections(delta_t)

    def update_infections(self, delta_t):
        """
        Computes new infections between lice and salmon. For each salmon within a
        fixed distance of a louse, we draw from an exp. dist. for a time. If the
        time is less than the timestep size, we infect.
        :param delta_t:
        :return:
        """
        # Get the un-attached lice
        unattached_lice = [louse for louse in self.lice if not louse.attached]

        # Calculate a louse-salmon distance matrix
        louse_positions = np.array([louse.position for louse in unattached_lice])
        louse_positions = louse_positions.reshape((len(unattached_lice), 1))
        salmon_positions = np.array([salmon.position for salmon in self.salmon])
        salmon_positions = salmon_positions.reshape((len(self.salmon), 1))

        dist_matrix = spatial.distance_matrix(louse_positions, salmon_positions)
        dist_matrix = np.abs(dist_matrix)  # We only care about the magnitude of distances

        # Find which salmon are within the distance threshold
        threshold = 3*self.salmon_sigma*np.sqrt(delta_t)
        threshold_matrix = dist_matrix < threshold

        # For each louse and each salmon within the threshold, check if we infect
        infections = 0

        for louse_idx in range(threshold_matrix.shape[0]):
            possible_targets = np.where(threshold_matrix[louse_idx, :])[0]

            # Sample from the exp. dist to check if we infect
            for target_salmon_idx in possible_targets:
                sampled_time = np.random.exponential(scale=1/self.infection_lambda)

                # If the salmon is infected, attach the louse to it
                #print(sampled_time, delta_t)
                if sampled_time <= delta_t:
                    #print('infected')
                    salmon = self.salmon[target_salmon_idx]
                    louse = unattached_lice[louse_idx]

                    louse.attach(salmon)
                    infections += 1
                    break
        #print(infections)

    def make_infection_plot(self, stage_to_plot=None, save_path=None):
        """
        Plots of histogram of louse positions for lice which are attached
        to a salmon. Can be restricted to certain stages of lice.

        :param stage_to_plot: The stage of lice we wish to count.
        :return:
        """
        infect_positions = []

        for salmon in self.salmon:
            if salmon.infected:
                for louse in salmon.lice:

                    # If we have specified a specific stage to plot, count only those lice
                    if stage_to_plot is not None:
                        if louse.stage == stage_to_plot:
                            infect_positions.append(salmon.position)
                    # Otherwise count all lice attached to a salmon
                    else:
                        infect_positions.append(salmon.position)

        plt.figure(figsize=(10, 10))

        if stage_to_plot is None:
            plt.title('Positions of attached lice: all stages')
        else:
            plt.title(f'Positions of attached lice: {stage_to_plot}')

        plt.xlabel('km downstream')
        plt.ylabel('Number of attached lice')

        # Bin the space into 1km blocks and create a normalized histogram of infections
        num_bins = int(self.end_x - self.start_x)
        plt.hist(infect_positions, histtype='step',
                 bins=np.linspace(self.start_x, self.end_x, num_bins), density=True)

        if save_path is not None:
            plt.savefig(save_path)

        plt.show()

    def make_abundance_plot(self, stage_to_plot, save_path=None):
        """
        Plots of abundance on each stage of lice per salmon. Abundance
        is the percentage of fish in a section of the river which

        :param stage_to_plot: The stage of lice we wish to count.
        :return:
        """
        # Bin the space into 1 km blocks
        num_bins = int(self.end_x - self.start_x)
        bins = np.linspace(self.start_x, self.end_x, num_bins, endpoint=True)
        infection_counts = np.zeros(num_bins)
        num_salmon = np.zeros(num_bins)

        for salmon in self.salmon:
            # Find the bin in which the count belongs
            for i in range(len(bins) - 1):
                if salmon.position >= bins[i] and salmon.position < bins[i+1]:
                    num_salmon[i] += 1

                    # See whether there are any attached lice of the stage
                    num_attached = 0

                    for louse in salmon.lice:
                        if louse.stage == stage_to_plot:
                            num_attached = 1
                            break

                    infection_counts[i] += num_attached
                    break

        # Normalize infection counts by number of salmon to get abundance
        abundances = np.divide(infection_counts, num_salmon)
        print(infection_counts, num_salmon)

        # Plot abundances
        plt.figure(figsize=(10, 10))
        plt.title(f'Relative abundance of {stage_to_plot}')
        plt.ylabel('Abundance')
        plt.xlabel('Position')
        plt.scatter(bins, abundances)
        plt.plot(bins, abundances)

        # Optionally save the figure
        if save_path is not None:
            plt.savefig(save_path)

        plt.show()

    def make_lice_counts_plot(self, stage_to_plot, save_path=None):
        """
        Plots of abundance on each stage of lice per salmon.

        :param stage_to_plot: The stage of lice we wish to count.
        :return:
        """
        infection_numbers = []

        for salmon in self.salmon:
            num_attached = 0

            for louse in salmon.lice:
                if louse.stage == stage_to_plot:
                    num_attached += 1

            infection_numbers.append(num_attached)

        # Plot abundances
        plt.figure(figsize=(10, 10))
        plt.hist(infection_numbers)

        if save_path is not None:
            plt.savefig(save_path)

        plt.show()

    def make_salmon_position_plot(self, save_path=None):
        """
        Plots the positions of each alive salmon.

        :return:
        """
        salmon_positions = []

        for salmon in self.salmon:
            salmon_positions.append(salmon.position)

        plt.figure(figsize=(10, 10))

        num_bins = int(self.end_x - self.start_x)
        plt.hist(salmon_positions, histtype='step', bins=np.linspace(self.start_x, self.end_x, num_bins))

        if save_path is not None:
            plt.savefig(save_path)

        plt.show()

    def make_louse_position_plot(self, save_path=None):
        """
        Plots the positions of all alive lice.

        :return:
        """
        louse_positions = []

        for louse in self.lice:
            louse_positions.append(louse.position)

        plt.figure(figsize=(10, 10))

        num_bins = int(self.end_x - self.start_x)
        plt.hist(louse_positions, histtype='step', bins=np.linspace(self.start_x, self.end_x, num_bins))

        if save_path is not None:
            plt.savefig(save_path)

        plt.show()

    def make_louse_age_plot(self, save_path=None):
        """
        Plots the histogram of lice ages.

        :return:
        """
        louse_ages = []

        for louse in self.lice:
            louse_ages.append(louse.age)

        plt.figure(figsize=(10, 10))
        plt.title('Histogram of Louse Ages')
        plt.hist(louse_ages, histtype='step')

        if save_path is not None:
            plt.savefig(save_path)

        plt.show()

def run_river(river, num_iters, delta_t):
    for i in tqdm(range(num_iters)):
        river.update(delta_t)

# river = River(salmon_velocity=1, salmon_sigma=0.5, louse_velocity=0.1, louse_sigma=0, salmon_spawn_rate=20,
#               louse_farm_rate=200, louse_ambient_rate=10, infection_lambda=2, end_x=20)
# river = River(salmon_velocity=0.3, salmon_sigma=0.3, louse_velocity=-0.1, louse_sigma=0.2589, salmon_spawn_rate=10,
#               louse_farm_rate=50, louse_ambient_rate=10, infection_lambda=0.00125, offspring_number=4, end_x=20)

river = River(salmon_velocity=0.3, salmon_sigma=0.3, louse_velocity=-0.1, louse_sigma=0.2589, salmon_spawn_rate=10,
              louse_farm_rate=50, louse_ambient_rate=10, infection_lambda=0.00125, offspring_number=4, end_x=20)

run_river(river, 1400, 0.1)

river.make_abundance_plot(stage_to_plot='copepodid')
#river.make_lice_counts_plot(stage_to_plot='copepodid')
# river.make_infection_plot(stage_to_plot='copepodid')
# river.make_salmon_position_plot()
river.make_louse_position_plot()
# river.make_louse_age_plot()
