def recalculate_order_cp2k(arguments):
    import argparse
    from types import SimpleNamespace
    
    import numpy as np
    import tomli
    
    from infretis.classes.engines.cp2k import read_xyz_file
    from infretis.classes.orderparameter import create_orderparameter
    
    parser = argparse.ArgumentParser(
        description="Recalculate the orderparameter from a .xyz file"
    )
    
    parser.add_argument("-xyz", help="The .xyz trajectory file")
    parser.add_argument("-box", help="box dimensions angstrom (e.g. 30 30 30)", type=float, nargs=3)
    parser.add_argument(
        "-toml", help="The .toml input file defining the orderparameter"
    )
    parser.add_argument(
        "-out",
        help="The output file. Default: order-rec.txt",
        default="order-rec.txt",
    )
    
    args = parser.parse_args(arguments)
    
    
    traj = read_xyz_file(args.xyz)
    box = np.array(args.box, dtype=float).flatten()
    with open(args.toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    
    orderparameter = create_orderparameter(toml_dict)
    # interfaces = toml_dict["simulation"]["interfaces"]
    
    with open(args.out, "w") as writefile:
        writefile.write("# step\torder\n")
        for i, frame in enumerate(traj):
            system = SimpleNamespace(
                pos=np.vstack((frame["x"], frame["y"], frame["z"])).T,
                box=box,
            )
            op = orderparameter.calculate(system)
            line = f"{i} " + " ".join([f"{opi}" for opi in op]) + "\n"
            writefile.write(line)
    
    print(f"\nAll done!\nOrderparameter values written to {args.out}")
