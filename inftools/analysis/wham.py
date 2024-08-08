from typing import Annotated
import typer

import os

import tomli

from inftools.analysis.Wham_Pcross import run_analysis


def wham(
    toml: Annotated[str, typer.Option("-toml", help="The infretis simulation toml file")] = "infretis.toml",
    data: Annotated[str, typer.Option("-data", help="The infretis data.txt file")] = "infretis_data.txt",
    nskip: Annotated[int, typer.Option("-nskip", help="Number of lines to skip in infretis_data.txt")] = 100,
    lamres: Annotated[float, typer.Option("-lamres", help="Resolution along the orderparameter, (intf1-intf0)/10)")] = None,
    nblock: Annotated[int, typer.Option("-nblock", case_sensitive=False, help="Minimal number of blocks in the block-error analysis")] = 5,
    fener: Annotated[bool, typer.Option("-fener", help="The infretis data.txt file")] = False,
    folder: Annotated[str, typer.Option("-folder", help="Output folder")] = "wham",
    ):
    """Run Titus0 wham script."""


    inps = {
        "toml": toml,
        "data": data,
        "nskip": nskip,
        "lamres": lamres,
        "nblock": nblock,
        "fener": fener,
        "folder": folder,
    }
    # load input:
    if os.path.isfile(inps["toml"]):
        with open(inps["toml"], mode="rb") as read:
            config = tomli.load(read)
    else:
        print("No toml file, exit.")
        return
    inps["intfs"] = config["simulation"]["interfaces"]

    if inps["lamres"] is None:
        inps["lamres"] = (inps["intfs"][1] - inps["intfs"][0]) / 10

    if inps["fener"]:
        inps["trajdir"] = config["simulation"]["load_dir"]

    run_analysis(inps)
