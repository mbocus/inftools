import os
from typing import Annotated

import matplotlib.pyplot as plt
import numpy as np
import tomli
import typer

from inftools.tistools.max_op import COLS


def plot_ens(
    toml: Annotated[str, typer.Option("-toml")] = "restart.toml",
    skip: Annotated[bool, typer.Option("-skip", help="skip initial load paths")] = False,
    save: Annotated[ str, typer.Option("-save", help="save with scienceplots.") ] = "no",
    pp: Annotated[
        bool, typer.Option("-pp", help="partial paths version")
    ] = False,
    cap: Annotated[
        int, typer.Option("-cap", help="max paths plotted per ens")
    ] = 100,
):
    """Plot sampled ensemble paths with interfaces"""
    if save != "no":
        plt.style.use("science")
        plt.figure(figsize=(14, 10))

    # Read toml info
    with open("restart.toml", "rb") as toml_file:
        toml = tomli.load(toml_file)
    intf = toml["simulation"]["interfaces"]
    datafile = toml["output"]["data_file"]
    load_dir = toml["simulation"]["load_dir"]

    plt.title("intfs: " + " ".join([str(i) for i in intf]))
    plt.axhline(intf[0], ls="--", color="k", alpha=0.5)
    plt.axhline(intf[-1], ls="--", color="k", alpha=0.5)

    cut = 5 if pp else 3

    acclen = 1
    for ens in list(range(len(intf))):
        cnt = 0
        if ens not in (0, len(intf) - 1):
            plt.axhline(intf[ens], color=f"{COLS[ens+1%len(COLS)]}", alpha=1.0)

        with open(datafile) as read:
            for line in read:
                if "#" in line:
                    continue
                rip = line.rstrip().split()
                frac = rip[cut : cut + len(intf)]
                if "-" in frac[ens]:
                    continue
                pnum = rip[0]
                pnum_f = f"{load_dir}/{pnum}/order.txt"
                if not os.path.isfile(pnum_f) or int(pnum) < len(intf):
                    continue
                data = np.loadtxt(pnum_f)
                plt.plot(data[:, 0] + acclen, data[:, 1], color=f"{COLS[ens%len(COLS)]}")
                acclen += len(data[:, 0])
                cnt += 1
                if cnt == cap:
                    break

    plt.xlabel(r"Time")
    plt.ylabel(r"Order Parameter")
    if save != "no":
        plt.savefig(save)
    else:
        plt.show()
