from typing import Annotated, Tuple
from enum import Enum
import typer

# Disable automatic underscore -> hyphen in CLI names
typer.main.get_command_name = lambda name: name

# for Choices CLI parameter
class Format(str, Enum):
    ONE = "trr"
    TWO = "xyz"
    THR = "traj"

def recalculate_order(
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    log: Annotated[str, typer.Option("-log")] = "sim.log",
    out: Annotated[bool, typer.Option("-out", help="the output of the analysis")] = False,
    format: Annotated[Format, typer.Option("-format", case_sensitive=False, help="the file format of the trajectory (.traj is ase format)")] = Format.ONE,
    box: Annotated[Tuple[int, int, int], typer.Option("-box", help="xyz only; box dimensions in angstrom (e.g. 30 30 30)")] = None,
    ):
    """
    Recalculate the orderparamter from a trajectory file.

    TODO:

    * We should be able to recalculate the OP from a path (i.e. use the
      traj.txt and the trajs in the accepted/ folder for some path)

    * Read velocity information, and set vel_rev in config
    """

    import os
    import numpy as np
    import tomli

    from infretis.classes.engines.cp2k import read_xyz_file
    from infretis.classes.engines.gromacs import read_trr_file
    from infretis.classes.orderparameter import create_orderparameter
    from infretis.classes.system import System

    if os.path.exists(out):
        raise FileExistsError(f"File {out} exists, aborting")

    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    orderparameter = create_orderparameter(toml_dict)
    # interfaces = toml_dict["simulation"]["interfaces"]

    if not traj.endswith(format):
        msg = f"[ WARNING ] {traj} may not be a {format} formatted"\
        " file. Use the -format flag to be sure you get some output"
        print(msg)
    if format == "xyz":
        if box is None:
            raise ValueError("Provide a box size for xyz files.")
        traj = read_xyz_file(traj)
        box = np.array(box, dtype=float).flatten()
    elif format == "traj":
        from ase.io.trajectory import Trajectory
        traj = Trajectory(traj)
    elif format == "trr":
        traj = read_trr_file(traj)


    with open(out, "w") as writefile:
        writefile.write("# step\torder\n")
        for i, frame in enumerate(traj):
            if format == "xyz":
                pos=np.vstack((frame["x"], frame["y"], frame["z"])).T
            elif format == "traj":
                pos = frame.positions
                box = np.diag(frame.cell.array)
            elif format == "trr":
                pos=frame[1]["x"]
                box=np.diag(frame[1]["box"])

            system = System()
            system.config = (traj, i)
            system.pos = pos
            system.box = box

            op = orderparameter.calculate(system)
            line = f"{i} " + " ".join([f"{opi}" for opi in op]) + "\n"
            writefile.write(line)

    print(f"[ INFO ] Orderparameter values written to {out}\n")

    return
