"""The functions to be used to run infretis via the terminal."""
import sys

from inftools.puckering_exercise.puckering import (
    check_indices,
    concatenate,
    generate_openff_topology,
    initial_path_from_iretis,
    initial_path_from_md,
    plot_order,
    recalculate_order,
)
from inftools.wham.wham import wham
from inftools.initial_paths import generate_zero_paths
from inftools.recalculate_order import recalculate_order_cp2k
from inftools.concatenate import concatenate_xyz
from inftools.get_interfaces import estimate_interfaces

# NOTE: when defining new functionality
# put the import statements in the function defenition
# as to avoid importing loads of libraries, which slows
# down the `inft` call from the command line

MAPPER = {
        "wham":wham,
        "check_indices":check_indices,
        "concatenate":concatenate,
        "concatenate_xyz":concatenate_xyz,
        "generate_openff_topology":generate_openff_topology,
        "initial_path_from_iretis":initial_path_from_iretis,
        "initial_path_from_md":initial_path_from_md,
        "plot_order":plot_order,
        "recalculate_order":recalculate_order,
        "recalculate_order_cp2k":recalculate_order_cp2k,
        "generate_zero_paths":generate_zero_paths,
        "est_intf": estimate_interfaces,
        }

def inftool():
    """Map an inftool command to a python function.

    Usage from the command line:
        inft `tool_name` `arguments`

    `tool_name`: the name of the tool you  want to use
    `arguments`: the arguments passed to the tool

    To get more information about each tool, pass `-h` after
    specifying the tool name, e.g. 'inft wham -h'.

    Example usage:
        inft wham -toml infretis.toml -data infretis_data.txt
        inft concatenate -h

    """
    if len(sys.argv) == 1 or sys.argv[1] in ["-h","help","--help"]:
        print(inftool.__doc__)
        print("Available inft commands:")
        for key in MAPPER.keys():
            print(f"\t{key}")
        return

    tool_name = sys.argv[1]
    arguments = sys.argv[2:]

    if tool_name not in list(MAPPER.keys()):
        msg = f"No tool named '{tool_name}', maybe you spelled it wrong?\n \
                \nFor usage information run\n\tinft -h"
        return msg

    tool = MAPPER[tool_name]
    # run the tool function
    tool(arguments)

def infretisinit():
    """To generate initial *toml template and other features."""
    return
