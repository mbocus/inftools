import argparse
from typing import Annotated
import typer


def plot_msg(
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    ):
    """Plot the order printed in the worker*/msg* files, so we
    can visualize the progress of the shooting moves.

    Require that the user be in the root sim folder."""
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from inftools.misc.infinit_helper import read_toml

    intfs = []
    # if toml, plot interfaces
    if os.path.isfile(toml):
        config = read_toml(toml)
        intfs = config["simulation"]["interfaces"]
        cap = config["simulation"]["tis_set"].get("interface_cap", None)
        plt.axhline(intfs[0], color='k')
        plt.axhline(intfs[-1], color='k')
        if cap is not None:
            plt.axhline(cap, color='r', ls="--")

    # iterate through worker folder
    for i in range(9999):
        wfolder = f"worker{i}"
        if not os.path.isdir(wfolder):
        	break
        idx = 0
        ens = -1

        # get op from msg files in sorted order
        files = [j for j in os.listdir(wfolder) if "msg" in j]
        files_idx = [int(j.split("_")[2]) for j in files]
        files = np.array(files)[np.argsort(files_idx)]

        for file in files:
            ens = int(file[4:7])
            msg = np.loadtxt(wfolder + "/" + file)
            try:
                plt.plot(msg[:, 0] + idx, msg[:, 1], color=f"C{i%8}")
                idx += len(msg[:, 0])
            except:
                print("cannot plot", wfolder + "/" + file, "perhaps no op printed yet")
        plt.scatter([msg[-1, 0] + idx -len(msg[:, 0])], [msg[-1, 1]], color=f"C{i%8}", edgecolor='k')

        # plt ensemble intfs
        if intfs and ens > -1:
            intf_idx = ens if ens == 0 else ens - 1
            plt.axhline(intfs[intf_idx], color=f"C{i%8}", label=f"{ens:03}-worker{i}")
            plt.legend()

    if i > 0:
        plt.show()
    else:
        print("Not in the root infretis simulation folder!")
