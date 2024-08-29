from typing import Annotated
import typer

import os

import tomli

from inftools.analysis.Wham_Pcross import run_analysis


def wham(
    toml: Annotated[str, typer.Option("-toml", help="The infretis .toml file")] = "infretis.toml",
    data: Annotated[str, typer.Option("-data", help="The infretis_data.txt file")] = "infretis_data.txt",
    nskip: Annotated[int, typer.Option("-nskip", help="Number of lines to skip in infretis_data.txt")] = 100,
    lamres: Annotated[float, typer.Option("-lamres", help="Resolution along the orderparameter, (intf1-intf0)/10)")] = None,
    nblock: Annotated[int, typer.Option("-nblock", case_sensitive=False, help="Minimal number of blocks in the block-error analysis")] = 5,
    fener: Annotated[bool, typer.Option("-fener", help="If set, calculate the conditional free energy. See Wham_")] = False,
    nbx: Annotated[int, typer.Option("-nbx", "Number of bins in x-direction when calculating the free-energy")] = 100,
    nby: Annotated[int, typer.Option("-nby", "Number of bins in x-direction when calculating the free-energy")] = 100,
    minx: Annotated[int, typer.Option("-minx", "Minimum orderparameter value in the x-direction when calculating FE")] = 0,
    maxx: Annotated[int, typer.Option("-maxx", "Maximum orderparameter value in the x-direction when calculating FE")] = 100,
    miny: Annotated[int, typer.Option("-miny", "Same as -minx but x-direction")] = 0,
    maxy: Annotated[int, typer.Option("-maxy", "Same as -maxx but x-direction")] = 360,
    xcol: Annotated[int, typer.Option("-xcol", "What column in order.txt to use as x-value when calculating FE")] = 1,
    ycol: Annotated[int, typer.Option("-ycol", "Same as -xcol but for y-value")] = 2,
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
        "histo_stuff":{
            "nbx":nbx, "minx":minx, "maxx":maxx, "xcol":xcol,
            "nby":nby, "miny":miny, "maxy":maxy, "ycol":ycol,
            }
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
