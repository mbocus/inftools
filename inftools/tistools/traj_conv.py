import argparse
from typing import Annotated
import typer

def trajtxt_conv(
    i: Annotated[str, typer.Option("-i")],
    o: Annotated[str, typer.Option("-o")],
    r: Annotated[int, typer.Option("-r", help="Subcycle ratio, subcycle xtc /subsycle trr")],
    ):
    """inft trjcat cannot be directly used when xtc is kept
    while delete_old=true during a gromacs infretis sim.

    Here we provide a quick traj.txt converter, that also consders
    if xtc frames are saved more frequently."""

    import numpy as np

    # Assert and data
    assert r >= 1
    data = np.loadtxt(i, dtype="str")
    data_out = []
    files = []

    # First get unique path files
    for file in data[:, 1]:
        if file not in files:
            files.append(file)

    # iterate over files
    for file_nr, file in enumerate(files):

        # get first, last and direction of frames in file
        file_idx_s = int(data[data[:, 1] == file][0][2])
        file_idx_e = int(data[data[:, 1] == file][-1][2])
        it = int(data[data[:, 1] == file][0][-1])

        # to account for overlap:
        # fill new data list
        extra_f, extra_b = 0, 0
        # usually first file trajectory contains all frames
        if file_nr > 0:
            # reverse,
            if it == -1:
                extra_b = 4
            elif it == 1 and file_idx_s == 1:
                extra_b = -4
                extra_f = 4
            # only one frame, traj usually contain at least two subcycles
            elif file_idx_s == file_idx_e:
                if it == 1:
                    extra_f = 4
                else:
                    extra_b = 4

        print(file)
        print(np.arange(file_idx_s*r + extra_b, file_idx_e*r + it + extra_f, it), len(np.arange(file_idx_s*r + extra_b, file_idx_e*r + it + extra_f, it)))
        print(file_idx_s, file_idx_e, it)
        print(file_idx_s*r + extra_b, file_idx_e*r + it + extra_f, extra_f, extra_b)
        print('')

        for idx in np.arange(file_idx_s*r + extra_b, file_idx_e*r + it + extra_f, it):
            data_out.append([file.replace(".trr", ".xtc"), idx, it])

    print("len", len(data_out))

    with open(o, "w") as write:
        write.write("# comment\n")
        write.write("# comment\n")
        for idx, line in enumerate(data_out):
            write.write(f"{idx:10.0f}    {line[0]}{line[1]:12.0f}{line[2]:7.0f}\n")

    print(f"trajtxt_conv Done! with file {o}")
