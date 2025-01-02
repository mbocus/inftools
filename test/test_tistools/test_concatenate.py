import os
import shutil

from inftools.tistools.concatenate import trjcat

def read_xyz(inp):
    cnt = 0
    with open(inp, 'r') as read:
        for line in read:
            if 'frame' in line:
                cnt += 1
    return cnt



def test_trjcat(tmp_path):
    """Check that we can modify the velocities with an engine,
    and that they are not equal to zero."""
    # folder we wil run from
    folder = tmp_path / "temp/"
    folder.mkdir()
    pathpath= os.path.dirname(__file__) + "/../simdata_dw/load/210"
    shutil.copytree(pathpath, str(folder) + "/path_xyz")
    os.chdir(str(folder)+"/path_xyz/accepted/")

    out = "merged.xyz"
    trjcat(out=out, traj="../traj.txt")
    cnt = read_xyz(out)

    # check if merged file exist and is of length 90
    assert os.path.isfile(out)
    assert cnt == 133

