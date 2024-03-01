from os import walk
import matplotlib.pyplot as plt
import numpy as np

def calc_rate(arguments):
    import argparse

    parser = argparse.ArgumentParser(
        description="Estimate interfaces from Pcross.txt"
    )
    parser.add_argument("-wham", help="The wham folder")
    parser.add_argument("-dt", type=float, help="timestep")
    parser.add_argument("-sc", type=float, help="subcycle")
    parser.add_argument("-unit", help="defaults: cp2k: 'fs', gromacs: 'ps'")
    parser.add_argument("-plot", help="specific value(s) vs runavg, in s-1.", nargs="+",)
    # parser.add_argument("-plot", type=float, action="store_true", help="Plot running rate with value")

    args = parser.parse_args(arguments)

    run_rate = np.loadtxt(args.wham + '/runav_rate.txt')
    run_rate[:, 3] = run_rate[:, 3]/(args.dt*args.sc)
    print(run_rate[-1, 3], f"{args.unit}^-1")
    if args.unit == 'fs':
        run_rate[:, 3] = run_rate[:, 3]*10**(15)
        print(run_rate[-1, 3], "s^-1")

    if args.plot:
        plt.plot(run_rate[:, 0], run_rate[:, 3])
        for i in args.plot:
            plt.axhline(float(i), color='k')
            minval = min(run_rate[-1, 3], float(i))
            maxval = max(run_rate[-1, 3], float(i))
            print('diff:', maxval/minval)
            # plt.plot([run_rate[0, 0], run_rate[-1, 0]], [i]
        plt.xlabel('cycles')
        plt.ylabel('s-1')
        plt.show()




