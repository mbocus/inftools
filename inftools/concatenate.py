def concatenate_xyz(arguments):
    import numpy as np
    import MDAnalysis as mda
    import os
    import argparse

    parser = argparse.ArgumentParser(
            description="Reverse and concatenate .xyz trajectories from an infretis simulation.")

    parser.add_argument("-out", help = "the outfile trajectory name")
    parser.add_argument("-traj", help="the traj.txt file. Trajectories should be in same folder")
    parser.add_argument("--selection", help="The selection, e.g. 'protein' or 'resname UNL' (default 'all')",
            default="all")

    args = parser.parse_args(arguments)

    traj_file_arr, index_arr = np.loadtxt(args.traj,usecols=[1,2],comments="#",dtype=str,unpack=True)
    index_arr = index_arr.astype(int)

    U = {}
    for traj_file in np.unique(traj_file_arr):
        if not os.path.exists(traj_file):
            raise FileNotFoundError(
                    f"\n No such file {traj_file}, maybe you forgot to 'cd accepted/'?")
        U[traj_file] = mda.Universe(traj_file)

    with mda.Writer(args.out,U[traj_file].select_atoms(args.selection).n_atoms) as wfile:
        for traj_file, index in zip(traj_file_arr,index_arr):
            u = U[traj_file]
            ag = u.select_atoms(args.selection)
            u.trajectory[index]
            wfile.write(ag.atoms)
    print("\nAll done!")
