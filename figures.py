import matplotlib.pyplot as plt

from agent_model import River
from tqdm import tqdm
import pandas as pd
import numpy as np
from scipy import interpolate


def run_river(river, num_iters, delta_t):
    """
    A helper function to iterate over an instance of a river

    :param river: The river to iterate over
    :param num_iters: The number of iterations
    :param delta_t: The timestep
    :return:
    """
    for i in tqdm(range(num_iters)):
        river.update(delta_t)

# Run the parameter sweep to find the optimal salmon velocity and lambda
def run_sweep_iter(salmon_velocity, lam, stage_to_plot='copepodid'):
    SPAWN_RATE = 20
    TIME_STEP = 0.1

    river = River(salmon_velocity=salmon_velocity, salmon_sigma=salmon_velocity, louse_velocity=-0.1,
                  louse_sigma=np.sqrt(0.2589), salmon_spawn_rate=SPAWN_RATE,
                  louse_farm_rate=5 * SPAWN_RATE, louse_ambient_rate=SPAWN_RATE,
                  infection_lambda=lam, offspring_number=17.57, end_x=20)

    # Compute how many iterations we need for salmon to make it to 40 km and run till then
    num_steps = int(40/(TIME_STEP*salmon_velocity))
    run_river(river, num_steps, TIME_STEP)

    # Read the true data and plot against the model
    if stage_to_plot == 'copepodid':
        true_df = pd.read_csv('data/copepodid_apr_1.csv', names=['x', 'y'])
        bins, abundances = river.make_abundance_plot(stage_to_plot='copepodid',
                                  save_path=f'sweep_out/{salmon_velocity}_velocity_{lam}_lambda_copepodid.pdf',
                                  show=False, true_data=true_df)
    elif stage_to_plot == 'chalimus':
        true_df = pd.read_csv('data/chalimus_apr_1.csv', names=['x', 'y'])
        bins, abundances = river.make_abundance_plot(stage_to_plot='chalimus',
                                                     save_path=f'sweep_out/{salmon_velocity}_velocity_{lam}_lambda_chalimus.pdf',
                                                     show=False, true_data=true_df)
    elif stage_to_plot == 'motiles':
        true_df = pd.read_csv('data/motile_apr_1.csv', names=['x', 'y'])
        bins, abundances = river.make_abundance_plot(stage_to_plot='motile',
                                                     save_path=f'sweep_out/{salmon_velocity}_velocity_{lam}_lambda_motile.pdf',
                                                     show=False, true_data=true_df)

    # Fit a smoothing spline
    spline = interpolate.CubicSpline(bins[~np.isnan(abundances)], abundances[~np.isnan(abundances)])

    # Evaluate the spline at the datapoint locations
    spline_eval = spline(true_df['x'].to_numpy())

    # Get distances to the true values
    dists = np.abs(spline_eval - true_df['y'].to_numpy())
    dist_avg = np.average(dists)

    return dist_avg

def make_abundances_panel(river):
    """
    Given a river, makes the relative abundances figure for all three
    stages of louse as side-by-side panels.

    :param river: The iterated river
    :return:
    """
    # Get the abundances for the three stages
    copepodid_bins, copepodid_abundances = river.make_abundance_plot(stage_to_plot='copepodid', show=False)
    chalimus_bins, chalimus_abundances = river.make_abundance_plot(stage_to_plot='chalimus', show=False)
    motile_bins, motile_abundances = river.make_abundance_plot(stage_to_plot='motile', show=False)

    # Fit a spline to the models
    copepodid_spline = interpolate.CubicSpline(copepodid_bins[~np.isnan(copepodid_abundances)],
                                               copepodid_abundances[~np.isnan(copepodid_abundances)])
    chalimus_spline = interpolate.CubicSpline(chalimus_bins[~np.isnan(chalimus_abundances)],
                                               chalimus_abundances[~np.isnan(chalimus_abundances)])
    motile_spline = interpolate.CubicSpline(motile_bins[~np.isnan(motile_abundances)],
                                               motile_abundances[~np.isnan(motile_abundances)])

    # Evaluate the splines at higher resolution
    new_bins = np.linspace(river.start_x, river.end_x, 200)
    copepodid_spline_eval = copepodid_spline(new_bins)
    chalimus_spline_eval = chalimus_spline(new_bins)
    motile_spline_eval = motile_spline(new_bins)

    # Get the ground truth data
    copepodid_gt = pd.read_csv('data/copepodid_apr_1.csv', names=['x', 'y'])
    chalimus_gt = pd.read_csv('data/chalimus_apr_1.csv', names=['x', 'y'])
    motile_gt = pd.read_csv('data/motile_apr_1.csv', names=['x', 'y'])

    plt.figure(figsize=(30, 10))

    plt.subplot(131)
    plt.title('Copepodid', fontsize=50)
    plt.xlim(river.start_x, river.end_x)
    plt.ylim(0, 1)
    plt.xticks(fontsize=30, rotation=45)
    plt.yticks(fontsize=30)
    plt.ylabel('Relative Abundance', fontsize=50)
    plt.plot(new_bins, copepodid_spline_eval, color='green')
    plt.plot(copepodid_gt['x'], copepodid_gt['y'], marker='o', color='red')

    plt.subplot(132)
    plt.title('Chalimus', fontsize=50)
    plt.xlim(river.start_x, river.end_x)
    plt.ylim(0, 1)
    plt.xticks(fontsize=30, rotation=45)
    plt.tick_params(
        axis='y',  # changes apply to the y-axis
        which='both',  # both major and minor ticks are affected
        left=False,  # ticks along the bottom edge are off
        labelleft=False)  # labels along the bottom edge are off
    plt.xlabel('Position in River (km)', fontsize=50)
    plt.plot(new_bins, chalimus_spline_eval, color='green')
    plt.plot(chalimus_gt['x'], chalimus_gt['y'], marker='o', color='red')

    plt.subplot(133)
    plt.title('Motile', fontsize=50)
    plt.xlim(river.start_x, river.end_x)
    plt.ylim(0, 1)
    plt.xticks(fontsize=30, rotation=45)
    plt.tick_params(
        axis='y',  # changes apply to the y-axis
        which='both',  # both major and minor ticks are affected
        left=False,  # ticks along the bottom edge are off
        labelleft=False)  # labels along the bottom edge are off

    plt.plot(new_bins, motile_spline_eval, label='Model Fit', color='green')
    plt.plot(motile_gt['x'], motile_gt['y'], marker='o', label='April 1st Dataset', color='red')
    plt.legend(fontsize=40)

    plt.tight_layout()
    plt.savefig('figures/relative_abundance_panel.pdf')

def make_planktonic_copepodid_comp(river_both, river_farm_only):
    """
    Makes a figure comparing the relative density of larval copepodids
    in a river which only has farm sources, as compared to one that has
    farm and ambient sources of lice.

    :param river_both: The iterated river with both sources.
    :param river_farm_only: The iterated river with only farm sources.
    :return:
    """
    # Get the planktonic louse positions for both rivers
    lice_both, bins_both = river_both.make_louse_position_plot(planktonic_only=True, show_plot=False)
    lice_farm, bins_farm = river_farm_only.make_louse_position_plot(planktonic_only=True, show_plot=False)

    # Get the midpoint of the bins to make plotting more convenient
    bins_both_mid = np.zeros(len(bins_both) - 1)
    bins_farm_mid = np.zeros(len(bins_farm) - 1)

    for i in range(len(bins_both) - 1):
        bins_both_mid[i] = (bins_both[i] + bins_both[i+1])/2
        bins_farm_mid[i] = (bins_farm[i] + bins_farm[i+1])/2

    bins_both = bins_both_mid
    bins_farm = bins_farm_mid

    # Fit splines to both
    both_spline = interpolate.CubicSpline(bins_both[~np.isnan(lice_both)],
                                               lice_both[~np.isnan(lice_both)])
    farm_spline = interpolate.CubicSpline(bins_farm[~np.isnan(lice_farm)],
                                          lice_farm[~np.isnan(lice_farm)])

    # Evaluate the splines
    new_bins_both = np.linspace(river_both.start_x, np.max(bins_both[~np.isnan(lice_both)]), 200)
    new_bins_farm = np.linspace(river_farm_only.start_x, np.max(bins_farm[~np.isnan(lice_farm)]), 200)
    both_spline_eval = both_spline(new_bins_both)
    farm_spline_eval = farm_spline(new_bins_farm)

    # Plot both results
    plt.figure(figsize=(10, 10))
    plt.ylabel('Planktonic Copepodid Density', fontsize=40)
    plt.xlabel('Position in River (km)', fontsize=40)
    plt.xlim(river_both.start_x, river_both.end_x)
    plt.ylim(0, 0.8)
    plt.xticks(fontsize=30, rotation=45)
    plt.yticks(fontsize=30)

    plt.plot(new_bins_both, both_spline_eval, color='blue', linewidth=4, label='Both Sources')
    plt.plot(new_bins_farm, farm_spline_eval, color='red', linewidth=4, label='Farm Only')

    plt.legend(fontsize=40)
    plt.tight_layout()
    plt.savefig('figures/planktonic_comparison.pdf')

def get_infection_pressure(river_both, river_ambient):
    """
    Computes the infection pressure in a river with both sources against that
    of a river with only ambient lice.

    :param river_both: A river with both sources of lice.
    :param river_ambient: A river with only ambient lice.
    :return: infection pressure (float)
    """
    # Find the number of infected salmon in the two rivers
    infected_both = 0
    infected_ambient = 0

    for salmon in river_both.salmon:
        if salmon.infected:
            infected_both += 1

    for salmon in river_ambient.salmon:
        if salmon.infected:
            infected_ambient += 1

    return infected_both/infected_ambient


def sweep_wrapper(args):
    """
    A multiprocessing wrapper for the run_sweep_iter function.

    :param args: An array of two values: velocity, lambda
    :return:
    """
    vel = args[0]
    lam = args[1]
    dist_avg, river = run_sweep_iter(vel, lam)
    return [vel, lam, dist_avg]

def run_sweep():
    """
    Runs the parameter sweep over a reasonable range of salmon velocities and lambdas.

    :return:
    """
    salmon_vel_range = np.linspace(0.45, 0.5, 1, endpoint=True)
    lambda_range = np.linspace(0.001, 0.0001, 1, endpoint=True)

    args_arr = []

    for vel in salmon_vel_range:
        for lam in lambda_range:
            args_arr.append([vel, lam])

    if __name__ == '__main__':
        with multiprocessing.Pool(5) as pool:
            results = list(pool.map(sweep_wrapper, args_arr))

        np.savetxt('sweep_out/results.txt', np.array(results))

# run_sweep()

# Make main figures
def make_figures():
    """
    Makes the figures featured in the main text of the report.

    :return:
    """
    BEST_VELOCITY = 0.2
    BEST_LAMBDA = 0.001
    SPAWN_RATE = 20
    TIME_STEP = 0.1
    RIVER_LENGTH = 40

    # Run the model
    river = River(salmon_velocity=BEST_VELOCITY, salmon_sigma=BEST_VELOCITY, louse_velocity=-0.1,
                  louse_sigma=np.sqrt(0.2589), salmon_spawn_rate=SPAWN_RATE,
                  louse_farm_rate=5 * SPAWN_RATE, louse_ambient_rate=SPAWN_RATE,
                  infection_lambda=BEST_LAMBDA, offspring_number=17.57, end_x=20)
    river_farm = River(salmon_velocity=BEST_VELOCITY, salmon_sigma=BEST_VELOCITY, louse_velocity=-0.1,
                  louse_sigma=np.sqrt(0.2589), salmon_spawn_rate=SPAWN_RATE,
                  louse_farm_rate=5 * SPAWN_RATE, louse_ambient_rate=0,
                  infection_lambda=BEST_LAMBDA, offspring_number=17.57, end_x=20)
    river_ambient = River(salmon_velocity=BEST_VELOCITY, salmon_sigma=BEST_VELOCITY, louse_velocity=-0.1,
                  louse_sigma=np.sqrt(0.2589), salmon_spawn_rate=SPAWN_RATE,
                  louse_farm_rate=0, louse_ambient_rate=SPAWN_RATE,
                  infection_lambda=BEST_LAMBDA, offspring_number=17.57, end_x=20)

    num_steps = int(RIVER_LENGTH/(TIME_STEP*BEST_VELOCITY))
    run_river(river, num_steps, TIME_STEP)
    run_river(river_farm, num_steps, TIME_STEP)
    run_river(river_ambient, num_steps, TIME_STEP)

    # Make the relative abundance panel
    make_abundances_panel(river)

    # Make the planktonic copepodid abundance plot (Krkosek fig. 3)
    make_planktonic_copepodid_comp(river, river_farm)

    # Get the infection pressure
    print('Infection pressure is:', get_infection_pressure(river, river_ambient))

make_figures()


