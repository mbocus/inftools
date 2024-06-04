def recalculate_order(arguments):
    """
    Recalculate the orderparamter from a trajectory file.

    TODO:
    * We should be able to recalculate the OP from a path (i.e. use the
      traj.txt and the trajs in the accepted/ folder for some path)

    * Read velocity information, and set vel_rev in config
    """

    import argparse

    import numpy as np
    import tomli

    from infretis.classes.engines.cp2k import read_xyz_file
    from infretis.classes.engines.gromacs import read_trr_file
    from infretis.classes.orderparameter import create_orderparameter
    from infretis.classes.system import System

    parser = argparse.ArgumentParser(
        description=recalculate_order.__doc__,
        formatter_class = argparse.RawTextHelpFormatter
    )

    parser.add_argument("-traj", help="The trajectory file")

    parser.add_argument("-box",
            help="xyz only; box dimensions in angstrom (e.g. 30 30 30)", type=float, nargs=3)

    parser.add_argument(
        "-toml", help="The .toml input file defining the orderparameter"
    )

    parser.add_argument("-format",
            help = "the file format of the trajectory (.traj is ase format)",
            choices = ["trr", "xyz", "traj"], default = "trr")

    parser.add_argument( "-out", help="The output file (order-rec.txt)", default="order-rec.txt")

    args = parser.parse_args(arguments)

    with open(args.toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    orderparameter = create_orderparameter(toml_dict)
    # interfaces = toml_dict["simulation"]["interfaces"]

    if args.format == "xyz":
        if args.box is None:
            raise ValueError("Provide a box size for xyz files.")
        traj = read_xyz_file(args.traj)
        box = np.array(args.box, dtype=float).flatten()
    elif args.format == "traj":
        from ase.io.trajectory import Trajectory
        traj = Trajectory(args.traj)
    elif args.format == "trr":
        traj = read_trr_file(args.traj)


    with open(args.out, "w") as writefile:
        writefile.write("# step\torder\n")
        for i, frame in enumerate(traj):
            if args.format == "xyz":
                pos=np.vstack((frame["x"], frame["y"], frame["z"])).T
            elif args.format == "traj":
                atoms = frame.positions
                box = np.diag(frame.box.cell.array)
            elif args.format == "trr":
                pos=frame[1]["x"]
                box=np.diag(frame[1]["box"])

            system = System()
            system.config = (args.traj, i)
            system.pos = pos
            system.box = box

            op = orderparameter.calculate(system)
            line = f"{i} " + " ".join([f"{opi}" for opi in op]) + "\n"
            writefile.write(line)

    print(f"\nAll done!\nOrderparameter values written to {args.out}")
