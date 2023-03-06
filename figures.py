from agent_model import River
from tqdm import tqdm
import pandas as pd

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
                  louse_sigma=0.2589, salmon_spawn_rate=SPAWN_RATE,
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
                  louse_sigma=0.2589, salmon_spawn_rate=SPAWN_RATE,
                  louse_farm_rate=5 * SPAWN_RATE, louse_ambient_rate=SPAWN_RATE,
                  infection_lambda=BEST_LAMBDA, offspring_number=17.57, end_x=20)

    num_steps = int(RIVER_LENGTH/(TIME_STEP*BEST_VELOCITY))
    run_river(river, num_steps, TIME_STEP)

    # Read the ground truth data files
    copepodid_df = pd.read_csv('data/copepodid_apr_1.csv', names=['x', 'y'])
    chalimus_df = pd.read_csv('data/chalimus_apr_1.csv', names=['x', 'y'])
    motile_df = pd.read_csv('data/motile_apr_1.csv', names=['x', 'y'])

    # Get the relative abundance figures
    bins, abundances = river.make_abundance_plot(stage_to_plot='copepodid',
                                                 save_path=f'figures/{BEST_VELOCITY}_velocity_{BEST_LAMBDA}_lambda_copepodid.pdf',
                                                 show=False, true_data=copepodid_df)
    bins, abundances = river.make_abundance_plot(stage_to_plot='chalimus',
                                                 save_path=f'figures/{BEST_VELOCITY}_velocity_{BEST_LAMBDA}_lambda_chalimus.pdf',
                                                 show=False, true_data=chalimus_df)
    bins, abundances = river.make_abundance_plot(stage_to_plot='motile',
                                                 save_path=f'figures/{BEST_VELOCITY}_velocity_{BEST_LAMBDA}_lambda_motile.pdf',
                                                 show=False, true_data=motile_df)

make_figures()


