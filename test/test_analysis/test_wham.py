import os
import shutil
import numpy as np

from inftools.analysis.wham import wham


def test_wham(tmp_path):
    """Check that we can modify the velocities with an engine,
    and that they are not equal to zero."""
    # folder we wil run from
    folder = tmp_path / "temp/"
    folder.mkdir()
    pathpath= os.path.dirname(__file__) + "/../simdata_dw/"
    shutil.copytree(pathpath, str(folder) + "/simdata")
    os.chdir(str(folder) + "/simdata")

    # run wham analysis script
    wham()

    data = np.loadtxt('wham/Pcross.txt')
    assert len(data) == 106
    assert data.shape == (106, 4)
    assert abs(6.958993086416705e-06 -data[-1][-1]) < 1000
