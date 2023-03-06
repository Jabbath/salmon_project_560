from agent_model import *

def salmon_test():
    test_salmon = Salmon(0, 5, 2)
    delta_t = 1

    for i in range(10):
        test_salmon.update_position(delta_t)
        print(test_salmon.position)

def louse_test():
    test_louse = Louse(0, .5, .2)
    delta_t = 1

    for i in range(10):
        test_louse.update_position(delta_t)
        print(test_louse.position)

def river_test():
    SPAWN_RATE = 10

    river = River(salmon_velocity=0.3, salmon_sigma=0.3, louse_velocity=-0.1,
                  louse_sigma=0.2589, salmon_spawn_rate=SPAWN_RATE,
                  louse_farm_rate=5*SPAWN_RATE, louse_ambient_rate=SPAWN_RATE,
                  infection_lambda=0.00125, offspring_number=17.57, end_x=20)
    run_river(river, 1400, 0.1)

    #river.make_abundance_plot(stage_to_plot='copepodid')
    #river.make_lice_counts_plot(stage_to_plot='copepodid')
    #river.make_infection_plot(stage_to_plot='chalimus')
    # river.make_salmon_position_plot()
    river.make_louse_position_plot()
    # river.make_louse_age_plot()