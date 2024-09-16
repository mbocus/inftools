import os
import pathlib
import typer

import numpy as np

from typing import Annotated
from infretis.classes.engines.factory import create_engines
from infretis.classes.orderparameter import create_orderparameters
from infretis.classes.path import Path, paste_paths
from infretis.classes.repex import REPEX_state
from infretis.classes.system import System
from infretis.setup import setup_config
from infretis.scheduler import scheduler

from inftools.analysis.infinit_helpers import *

PHRASES = [
    ["Infinitely fast from A to B","with ∞RETIS"],
    ["Infinit mode:", "Engaging endless loops with ∞RETIS",],
    ["Installing infinite improbability drivers ...", "∞RETIS is allready installed!",],
    ["Solving the time warp", "in digital chemical discoveries",],
    ["Performing infinite swaps", "no molecules are harmed in this exchange!",],
    ["Asynchronous swaps engaged", "because who has time for synchronization?",],
    ["Pushing for transitions", "molecules, please keep your seatbelts fastened!",],
    ["Propagating through time", " please hold on to your atoms!",],
    ["Propagating backwards in time", "because forward is too mainstream!",],
    ["Fluxing through rare events",  "with ∞RETIS",],
    ["Taking molecular strolls in parallel", "with ∞RETIS",],
    ["Shooting through the void", "with ∞RETIS",],
    ["Infinit /ˈɪnfɪnət/ (adjective)", "limitless or endless in space, extent, or size"]
]


def generate_zero_paths(
    conf: Annotated[str, typer.Option("-conf", help="The name (not the path) of the initial configuration to propagate from. inftools will look in the input folder specified in the .toml file.")],
    toml: Annotated[str, typer.Option("-toml",)] = "infretis.toml",
    config : Annotated[str, typer.Option(hidden = True)] = None,
    ):
    """Generate initial paths for the [0-] and [0+]

    ensembles by propagating a single configuration forward
    and backward in time until it crosses the lambda0 interface.
    These can be used to e.g. push the system up the barrier using
    multiple infRETIS simulations."""

    # make a directory we work from
    tmp_dir = pathlib.Path("temporary_load/")
    tmp_dir.mkdir(exist_ok = False)
    load_dir = pathlib.Path("load/")
    load_dir.mkdir(exist_ok = False)

    initial_configuration = conf

    # maximal length of initial paths

    # infretis parameters
    if config == None:
        config = setup_config(toml)
    maxlen = config["simulation"]["tis_set"]["maxlength"]
    state = REPEX_state(config, minus=True)

    # setup ensembles
    state.initiate_ensembles()
    state.engines = create_engines(config)
    create_orderparameters(state.engines, config)

    # initial configuration to start from
    system0 = System()
    engine = state.engines[0]
    engine.exe_dir = str(tmp_dir.resolve())
    if "dask" in config.keys():
        wmdrun = config["dask"]["wmdrun"][0]
    else:
        wmdrun = config["runner"]["wmdrun"][0]
    engine.set_mdrun(
        {"wmdrun": wmdrun, "exe_dir": engine.exe_dir}
    )
    system0.set_pos((os.path.abspath(initial_configuration), 0))
    engine.rgen = np.random.default_rng()
    engine.modify_velocities(system0, config["simulation"]["tis_set"])

    # empty paths we will fill forwards in time in [0-] and [0+]
    path0 = Path(maxlen=maxlen)
    path1 = Path(maxlen=maxlen)

    # propagate forwards from the intiial configuration
    # note that one of these does not get integrated because
    # the initial phasepoint is either below or above interface 0
    print("Propagating in ensemble [0-]")
    status0, message0 = engine.propagate(path0, state.ensembles[0], system0)
    system0.set_pos((os.path.abspath(initial_configuration), 0))
    print("Propagating in ensemble [0+]")
    status1, message1 = engine.propagate(path1, state.ensembles[1], system0)

    # we did only one integration step in ensemble 0 because
    # we started above interface 0
    if path0.length == 1:
        print("Re-propagating [0-] since we started above lambda0")
        system0.set_pos((engine.dump_config(path1.phasepoints[-1].config), 0))
        path0 = Path(maxlen=maxlen)
        status0, message0 = engine.propagate(path0, state.ensembles[0], system0)

    # or we did only one integration step in ensemble 1 because
    # we started below interface 0
    elif path1.length == 1:
        print("Re-propagating [0+] since we started below lambda0")
        system0.set_pos((engine.dump_config(path0.phasepoints[-1].config), 0))
        path1 = Path(maxlen=maxlen)
        status1, message1 = engine.propagate(path1, state.ensembles[1], system0)

    else:
        raise ValueError("Something fishy!\
                Path lengths in one of the ensembles != 1")

    # backward paths
    path0r = Path(maxlen=maxlen)
    path1r = Path(maxlen=maxlen)

    print("Propagating [0-] in reverse")
    status0, message0 = engine.propagate(
        path0r, state.ensembles[0], path0.phasepoints[0], reverse=True
    )

    print("Propagating [0+] in reverse")
    status1, message1 = engine.propagate(
        path1r, state.ensembles[1], path1.phasepoints[0], reverse=True
    )

    print("Done! Making load/ dir")
    # make load directories
    pathsf = [path0, path1]
    pathsr = [path0r, path1r]
    for i in range(2):
        dirname = load_dir / str(i)
        accepted = dirname / "accepted"
        orderfile = dirname / "order.txt"
        trajtxtfile = dirname / "traj.txt"
        print(f"Making folder: {str(dirname)}")
        dirname.mkdir()
        print(f"Making folder: {str(accepted)}")
        accepted.mkdir()
        # combine forward and backward path
        path = paste_paths(pathsr[i], pathsf[i])
        # save order paramter
        order = [pp.order[0] for pp in path.phasepoints]
        order = np.vstack((np.arange(len(order)), np.array(order))).T
        np.savetxt(str(orderfile), order, fmt=["%d", "%12.6f"])
        N = len(order)
        # save traj.txt
        np.savetxt(
            str(trajtxtfile),
            np.c_[
                [str(i) for i in range(N)],
                [pp.config[0].split("/")[-1] for pp in path.phasepoints],
                [pp.config[1] for pp in path.phasepoints],
                [-1 if pp.vel_rev else 1 for pp in path.phasepoints],
            ],
            header=f"{'time':>10} {'trajfile':>15} {'index':>10} {'vel':>5}",
            fmt=["%10s", "%15s", "%10s", "%5s"],
        )
        # copy paths
        for trajfile in np.unique(
            [pp.config[0].split("/")[-1] for pp in path.phasepoints]
        ):
            traj_path = tmp_dir / trajfile
            traj_path.rename(accepted / trajfile)


def infinit(
    toml: Annotated[str, typer.Option("-toml", help="Path to .toml")] = "infretis.toml",
    log: Annotated[str, typer.Option("-log", help="File for logging output")] = "infinit.log",
    ):
    """The infretis initial path generator.

    Based on the YouTube series:
    https://www.youtube.com/watch?v=mW9tC2A7COs&list=PL5dSi5eZMe1iN_Uz8pTph6i8AGXhVUZIj&index=24
    """

    # Lecture 04:
    # define the grid spacing for lambda values. All interface positions are
    # are rounded off to this vlue
    lamres = 0.01 # e.g. second digit if lamres = 0.01

    # skip this fraction of initial paths for analysis (when estimating new intf?)
    initskip = 0.1 # skip 10% of initial paths

    # compute efficiency measure for present set of interfaces and optimal set
    # (estimated from WHAM). If efficiency is worse than some factor, we update
    tolerance = 1.2 # if worse than 20% compared to optimal efficiency

    # estimated lower bound for local crossing probability
    pL = 0.3

    # njumps Lp/Ls where Lp is average len full path, and Ls average len sub traj
    # lambda_cap to avoid A -> B trajs. E.g. 10% B->B paths in wf. Alternatively,
    # palce lambda_cap in half between lambdaN and lambda(N-1)

    # Lecture 05: set lambda0 lambdaN


    # TODO: restart, log file instead of print
    log = LightLogger(log)


    # we need among others parameters set in [infinit]
    config = setup_config(toml)
    # get the infinit settings from 'config' and set default parameters
    iset = set_default_infinit(config)
    cstep = iset["cstep"]

    if iset["cstep"]  == 0:
        log.log("Generating zero paths ...")
        init_conf = pl.Path(iset["initial_conf"]).resolve()
        generate_zero_paths(str(init_conf), toml = toml, config = config)
        log.log("Done with zero paths!\n")
    log.log('Running infretis initialization "infinit" ...')
    print_logo()
    for iretis_steps in iset["steps_per_iter"][cstep:]:
        log.log(f"Step {iset['cstep']}:")
        #print_logo(step = iset["cstep"])
        #time.sleep(3)
        run_infretis(config, iretis_steps)
        iset["cstep"] += cstep + 1
