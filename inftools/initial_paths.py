import argparse
import os
import pathlib
import shutil

import numpy as np

from infretis.classes.engines.factory import create_engines
from infretis.classes.orderparameter import create_orderparameters
from infretis.classes.path import Path, paste_paths
from infretis.classes.repex import REPEX_state
from infretis.classes.system import System
from infretis.setup import setup_config


def generate_zero_paths(arguments):
    parser = argparse.ArgumentParser(
        description="Generate initial paths for the [0-] and [0+] \
                ensembles by propagating a single configuration forward \
                and backward in time until it crosses the lambda0 interface.\
                These can be used to e.g. push the system up the barrier using \
                multiple infRETIS simulations."
    )

    parser.add_argument("-toml", help = "The .toml file")
    parser.add_argument("-conf",
        help = "The initial configuration to propagate from")
    parser.add_argument("-maxlen", help = "The maximum allowed path length", type = int)

    args = parser.parse_args(arguments)

    # make a directory we work from
    tmp_dir = pathlib.Path("temporary_load/")
    tmp_dir.mkdir(exist_ok = False)
    load_dir = pathlib.Path("load/")
    load_dir.mkdir(exist_ok = False)

    initial_configuration = args.conf

    # maximal length of initial paths
    maxlen = args.maxlen
    
    # infretis parameters
    config = setup_config(args.toml)
    config["output"]["data_dir"] = str(tmp_dir)
    state = REPEX_state(config, minus=True)
    
    # setup ensembles
    state.initiate_ensembles()
    state.engines = create_engines(config)
    create_orderparameters(state.engines, config)
    
    # initial configuration to start from
    system0 = System()
    engine = state.engines[config["engine"]["engine"]]
    engine.exe_dir = str(tmp_dir.resolve())
    print(engine.exe_dir)
    engine.set_mdrun(
        {"wmdrun": config["dask"]["wmdrun"][0], "exe_dir": engine.exe_dir}
    )
    system0.set_pos((os.path.join(engine.input_path, initial_configuration), 0))
    
    # empty paths we will fill forwards in time in [0-] and [0+]
    path0 = Path(maxlen=maxlen)
    path1 = Path(maxlen=maxlen)
    
    # propagate forwards from the intiial configuration
    # note that one of these does not get integrated because
    # the initial phasepoint is either below or above interface 0
    print("Propagating in ensemble [0-]")
    status0, message0 = engine.propagate(path0, state.ensembles[0], system0)
    system0.set_pos((os.path.join(engine.input_path, initial_configuration), 0))
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
            shutil.copy(tmp_load / trajfile, accepted)
