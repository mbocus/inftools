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


def inftool():
    """Map an inftool command to a python function.

    Usage from the command line
        inft `tool_name` `arguments`

    `tool_name` is a the name of the tool you 
    want to use.
    `arguments` the arguments passed to the tool.

    To get more information about each tool, pass `-h` after
    specifying the tool name.
    """

    tool_name = sys.argv[1]
    arguments = sys.argv[1:]
    
    # NOTE: when defining new functionality
    # put the import statements in the function defenition
    # as to avoid importing loads of libraries, which slows
    # down the `inft` call from the command line
    mapper = {
            "wham":wham,
            "check_indices":check_indices,
            "concatenate":concatenate,
            "generate_openff_topology":generate_openff_topology,
            "initial_path_from_iretis":initial_path_from_iretis,
            "initial_path_from_md":initial_path_from_md,
            "plot_order":plot_order,
            "recalculate_order":recalculate_order,
            }

    if sys.argv[1] in ["-h","help","--help"]:
        print(inftool.__doc__)
        print("Available inft commands:")
        for key in mapper.keys():
            print(f"\t{key}")
        return

    tool = mapper[tool_name]
    # run the tool function
    tool(arguments)

def infretisinit():
    """To generate initial *toml template and other features."""
    return
