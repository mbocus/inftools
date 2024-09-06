import struct

import numpy as np

from infretis.classes.engines.gromacs import (
read_trr_data,
skip_trr_data,
read_trr_header
)

def read_trr_file(filename: str, read_data: bool = True):
    """Yield frames from a TRR file."""
    with open(filename, "rb") as infile:
        while True:
            try:
                header, _ = read_trr_header(infile)
                if read_data:
                    data = read_trr_data(infile, header)
                else:
                    skip_trr_data(infile, header)
                    data = None
                yield header, data
            except EOFError:
                return None, None
            except struct.error:
                logger.warning(
                    "Could not read a frame from the TRR file. Aborting!"
                )
                return None, None
