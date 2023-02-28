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
#salmon_test()
#louse_test()