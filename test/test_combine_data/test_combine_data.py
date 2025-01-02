"""Test methods for combine_data."""
import os
from distutils.dir_util import copy_tree
from pathlib import PosixPath

import numpy as np
import pytest
import tomli

from inftools.analysis.wham import wham
from inftools.tistools.combine_results import combine_data


@pytest.mark.heavy
def test_infinit_1(tmp_path: PosixPath) -> None:
    """Test infinit from phase point"""
    folder = tmp_path / "temp"
    folder.mkdir()
    basepath = PosixPath(__file__).parent
    data_dir = (basepath / "data").resolve()
    copy_tree(str(data_dir), str(folder))
    os.chdir(folder)

    tomls=["sim1.toml", "sim2.toml", "sim3.toml"]
    combine_data(
        tomls=tomls,
        datas=["sim1.txt.gz", "sim2.txt.gz", "sim3.txt.gz"],
        skip=100,
        scramble=True,
        out="combo"
    )

    # Get the combo interfaces.
    with open("combo.toml", "rb") as rfile:
        conf = tomli.load(rfile)
        combo_intf = conf["simulation"]["interfaces"]

    intfs = []
    for toml in tomls:
        with open(toml, "rb") as rfile:
            conf = tomli.load(rfile)
            intf = conf["simulation"]["interfaces"]
            assert len(combo_intf) >= len(intf)
            intfs += intf
    assert set(intfs) == set(combo_intf)

    # Assert that files exist
    assert os.path.isfile("combo.toml")
    assert os.path.isfile("combo.txt")

    # Run wham
    wham(toml="combo.toml", data="combo.txt")

    # Assert that important files are created
    assert os.path.isfile("wham/runav_rate.txt")
    assert os.path.isfile("wham/runav_flux.txt")
    assert os.path.isfile("wham/runav_Pcross.txt")

    # Test that we get good results
    data = np.loadtxt("wham/runav_rate.txt")
    calc_rate = data[-1, -1]/0.025
    kramer = 2.58 * 10 **(-7) # from wf paper
    assert 100*abs(kramer-calc_rate)/kramer < 7 # %
