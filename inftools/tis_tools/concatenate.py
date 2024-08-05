from typing import Annotated
import typer

def concatenate_xyz(
    out: Annotated[str, typer.Option("-out", help="the outfile trajectory name")],
    traj: Annotated[str, typer.Option("-traj", help="the traj.txt file. Trajectories should be in same folder")],
    selection: Annotated[str, typer.Option("-selection", help="The selection, e.g. 'protein' or 'resname UNL'")]="all",
    engine: Annotated[str, typer.Option("-engine", help="The selection, e.g. 'protein' or 'resname UNL'")]="cp2k",
        ):
    "Reverse and concatenate .xyz trajectories from an infretis simulation."
    import numpy as np
    import MDAnalysis as mda
    import os
    import argparse
    from ase.io.trajectory import Trajectory
    from ase.io import write

    traj_file_arr, index_arr = np.loadtxt(traj,usecols=[1,2],comments="#",dtype=str,unpack=True)
    index_arr = index_arr.astype(int)

    # read the trajectories
    U = {}
    for traj_file in np.unique(traj_file_arr):
        if not os.path.exists(traj_file):
            raise FileNotFoundError(
                    f"\n No such file {traj_file}, maybe you forgot to 'cd accepted/'?")
        if engine == "cp2k":
            U[traj_file] = mda.Universe(traj_file)
            n_atoms = U[traj_file].select_atoms(selection).n_atoms
        elif engine == "ase":
            if selection != "all":
                raise NotImplementedError("Only selection all is implemented for ase egnine.")
            U[traj_file] = Trajectory(traj_file)
            n_atoms = U[traj_file][0].positions.shape[0]

    # write the concatenated trajectory
    if engine == "cp2k":
        with mda.Writer(out, n_atoms) as wfile:
            for traj_file, index in zip(traj_file_arr,index_arr):
                u = U[traj_file]
                ag = u.select_atoms(selection)
                u.trajectory[index]
                wfile.write(ag.atoms)
    elif engine == "ase":
        for traj_file, index in zip(traj_file_arr, index_arr):
            # traj object
            u = U[traj_file]
            atoms = u[index]
            write(out, atoms, append = True)

    print("\nAll done!")

