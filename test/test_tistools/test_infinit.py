"""Test methods for doing TIS."""
import difflib
import filecmp
import os
import shutil
from pathlib import PosixPath
from subprocess import STDOUT, check_output

import pytest
import tomli
import tomli_w


@pytest.mark.heavy
def test_infinit_1(tmp_path: PosixPath) -> None:
    """Test infinit from phase point"""
    folder = tmp_path / "temp"
    folder.mkdir()
    basepath = PosixPath(__file__).parent
    load_dir = (
        basepath / "../../examples/turtlemd/double_well/load_copy"
    ).resolve()
    toml_dir = basepath / "data/infretis.toml"
    conf_dir = basepath / "data/initial.xyz"
    shutil.copy(str(toml_dir), str(folder))
    shutil.copy(str(conf_dir), str(folder))
    os.chdir(folder)

    print('banda', folder)
    # assert False
    success = os.system("inft infinit")

    # Check if files exist and that p is same or growing
    tomls = [f"infretis_{i}.toml" for i in range(2, 11)] + ["infretis.toml"]
    assert os.path.isfile(tomls[0])
    config0 = tomli.load(tomls[0])
    assert "infinit" in config0
    assert "p" in config0["infinit"]
    for item in items[1:]:
        assert os.path.isfile(item)
        config = tomli.load(item)
        assert "infinit" in config
        assert "p" in config["infinit"]
        assert len(config["infinit"]["p"]) >= len(config0["infinit"]["p"])


@pytest.mark.heavy
def test_infinit_2(tmp_path: PosixPath) -> None:
    """Test infinit from existing sim"""
    folder = tmp_path / "temp"
    folder.mkdir()
    basepath = PosixPath(__file__).parent
    load_dir = (
        basepath / "../../examples/turtlemd/double_well/load_copy"
    ).resolve()
    toml_dir = basepath / "data/infretis.toml"
    conf_dir = basepath / "data/initial.xyz"
    shutil.copy(str(toml_dir), str(folder))
    shutil.copy(str(conf_dir), str(folder))
    os.chdir(folder)

    print('banda', folder)
    # assert False
    success = os.system("inft infinit")

    # Check if files exist and that p is same or growing
    tomls = [f"infretis_{i}.toml" for i in range(2, 11)] + ["infretis.toml"]
    assert os.path.isfile(tomls[0])
    config0 = tomli.load(tomls[0])
    assert "infinit" in config0
    assert "p" in config0["infinit"]
    for item in items[1:]:
        assert os.path.isfile(item)
        config = tomli.load(item)
        assert "infinit" in config
        assert "p" in config["infinit"]
        assert len(config["infinit"]["p"]) >= len(config0["infinit"]["p"])
