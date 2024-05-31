def calc_flow(arguments):
    """
    Keep track of which parent paths are in which ensemble after each MC move.

    We start the simulation with a valid path in each ensemble, and we
    would like to see how the initial paths move through the ensembles.

    For example, when starting out, path nr. 4 is initially in ensemble 004.
    However, it can be swapped with other ensembles, and as such move upward
    all the way til it become a reactive path, or it may move down all the way
    til it becomes a path in [0-].

    The output of this function returns a text file where the column i specifies
    which ensemble the parent path i is currently in.
    """
    import numpy as np
    import tomli
    import argparse
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(
        description=calc_flow.__doc__)

    parser.add_argument("-toml", default = "infretis.toml",
            help = ".toml file to read nr. of ensembles")
    parser.add_argument("-log", default="sim.log",
            help = "the .log file to read path numbers")
    parser.add_argument("-out", default=False,
            help = "the output of the analysis")
    parser.add_argument("-plot", nargs='*',
            help = "if integers, plot the flow for those paths, if 'all' plot flow for all paths")

    args = parser.parse_args(arguments)

    with open(args.toml, "rb") as toml_file:
        config = tomli.load(toml_file)

    n_ensembles = len(config["simulation"]["interfaces"])

    # map a path number from sim.log to its initial parent path
    path_map = {}
    # keep track of which parent path is in which ensemble
    flow_map = [[] for i in range(n_ensembles)]
    # if true we initialize by reading the first matrix in the sim.log
    # to get the parent path numbers
    initialize = True
    # we found the beginning of the first matrix in sim.log
    found_start_init = False
    # number of paths read from the first matrix in sim.log when initializing
    n_paths_read_init = 0

    with open(args.log, "r") as rfile:
        for i,line in enumerate(rfile):
            #print(i)
            if "->" in line and not initialize:
                line = line.split()
                # this is a zero swap
                if len(line)==15:
                    ens0, ens1 = 0, 1
                    pathnr0_old, pathnr1_old = line[-5:-3]
                    pathnr1_new, pathnr0_new = line[-2:]

                    if pathnr0_new != pathnr0_old:
                        path_map[pathnr0_new] = path_map[pathnr0_old]
                        del path_map[pathnr0_old]

                        path_map[pathnr1_new] = path_map[pathnr1_old]
                        del path_map[pathnr1_old]

                    pathnr_ens0 = path_map[pathnr0_new]
                    pathnr_ens1 = path_map[pathnr1_new]
                    flow_map[pathnr_ens0].append(1)
                    flow_map[pathnr_ens1].append(0)
                    for j in range(n_ensembles):
                        if j != pathnr_ens0 and j!= pathnr_ens1:
                            flow_map[j].append(flow_map[j][-1])

                elif len(line)==11:
                    pathnri_old = line[-3]
                    pathnri_new = line[-1]
                    ensi = int(line[-6])

                    if pathnri_new != pathnri_old:
                        path_map[pathnri_new] = path_map[pathnri_old]
                        del path_map[pathnri_old]

                    pathnr_ensi = path_map[pathnri_new]

                    for j in range(n_ensembles):
                        if flow_map[j][-1] == ensi and j != pathnr_ensi:
                            flow_map[j].append(flow_map[pathnr_ensi][-1])
                        elif j != pathnr_ensi:
                            flow_map[j].append(flow_map[j][-1])

                    flow_map[pathnr_ensi].append(ensi)

                else:
                    ValueError("Wrong length in line.split()")
                #print("=="*10)
                #print(line)
                #print(path_map)
                #print(flow_map)
                #print("=="*10)

            elif initialize and "  -- |" in line:
                found_start_init = True

            elif initialize and found_start_init:
                line = line.split()
                pathnri = str(int(line[1][1:]))
                path_map[pathnri] = n_paths_read_init
                flow_map[n_paths_read_init].append(n_paths_read_init)
                n_paths_read_init += 1
                if n_paths_read_init == n_ensembles:
                    initialize = False
                    print(f"[INFO] Initialization done {path_map}")

            else:
                ValueError("Something wrong, but not sure what")

        # write to file
        if args.out:
            np.savetxt(args.out, np.array(flow_map).T, fmt="%4d")


        # plotting stuff
        if args.plot is not None:
            if args.plot == ['all']:
                plot_range = range(n_ensembles)
            else:
                plot_range = []
                for i in args.plot:
                    try:
                        plot_range.append(int(i))
                    except ValueError:
                        print(f"Invalid axis value: {args.plot}")
            for j in plot_range:
                plt.plot(flow_map[j], label = f"path {j}", marker = "o", markersize = 2)
            plt.legend()
            plt.show()
