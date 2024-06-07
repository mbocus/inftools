def concatenate_xyz(arguments):
    import numpy as np
    import MDAnalysis as mda
    import os
    import argparse
    from ase.io.trajectory import Trajectory
    from ase.io import write

    parser = argparse.ArgumentParser(
            description="Reverse and concatenate .xyz trajectories from an infretis simulation.")

    parser.add_argument("-out", help = "the outfile trajectory name")
    parser.add_argument("-traj", help="the traj.txt file. Trajectories should be in same folder")
    parser.add_argument("--selection", help="The selection, e.g. 'protein' or 'resname UNL' (default 'all')",
            default="all")
    parser.add_argument("-engine", help = "The engine used to create the trajectories", default = "cp2k")

    args = parser.parse_args(arguments)


    traj_file_arr, index_arr = np.loadtxt(args.traj,usecols=[1,2],comments="#",dtype=str,unpack=True)
    index_arr = index_arr.astype(int)

    # read the trajectories
    U = {}
    for traj_file in np.unique(traj_file_arr):
        if not os.path.exists(traj_file):
            raise FileNotFoundError(
                    f"\n No such file {traj_file}, maybe you forgot to 'cd accepted/'?")
        if args.engine == "cp2k":
            U[traj_file] = mda.Universe(traj_file)
            n_atoms = U[traj_file].select_atoms(args.selection).n_atoms
        elif args.engine == "ase":
            if args.selection != "all":
                raise NotImplementedError("Only selection all is implemented for ase egnine.")
            U[traj_file] = Trajectory(traj_file)
            n_atoms = U[traj_file][0].positions.shape[0]

    # write the concatenated trajectory
    if args.engine == "cp2k":
        with mda.Writer(args.out, n_atoms) as wfile:
            for traj_file, index in zip(traj_file_arr,index_arr):
                u = U[traj_file]
                ag = u.select_atoms(args.selection)
                u.trajectory[index]
                wfile.write(ag.atoms)
    elif args.engine == "ase":
        for traj_file, index in zip(traj_file_arr, index_arr):
            # traj object
            u = U[traj_file]
            atoms = u[index]
            write(args.out, atoms, append = True)

    print("\nAll done!")

