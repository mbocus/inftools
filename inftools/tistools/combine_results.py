from typing import Annotated

import numpy as np
import tomli
import tomli_w
import typer

from inftools.misc.data_helper import data_reader


def combine_data(
    tomls: Annotated[list[str], typer.Option("-tomls", help="tomls for all \
simulations")],
    datas: Annotated[list[str], typer.Option("-datas", help="data files for \
all simulations")],
    skip: Annotated[int, typer.Option("-skip", help="skip initial lines for \
all simulations")] = 100,
    scramble: Annotated[bool, typer.Option("-scramble", help="If output combo \
data lines are scrambled.")] = True,
    out: Annotated[str, typer.Option("-out", help="name for output .txt/toml \
file.")] = "combo"
):
    """Combine different infretis simulations.

    A relatively convoluted task given that we have to merge
    simulation ensembles together from small to big (if I
    assume correctly).

    Require input of their individual toml and data files,
    in their respective order, e.g.

    -tomls sim1 -tomls sim2 -datas data1 -datas data2.

    # NB: If output data is not scrambled, then one
    column will remain as "----" state for one whole sim.
    On the other side, in this implementation all data must
    then be stored in ram.
    """
    # do some initial checks
    assert len(set(tomls)) == len(tomls)
    assert len(set(datas)) == len(datas)

    # initialize some variables
    sims = {}
    intfa, intfb, intfs = [], [], []

    # loop through individual simulations
    for idx, (toml, data) in enumerate(zip(tomls, datas)):
        sims[idx] = {"toml": toml, "data": data}
        with open(toml, "rb") as rfile:
            sims[idx]["conf"] = tomli.load(rfile)
            sims[idx]["intf"] = sims[idx]["conf"]["simulation"]["interfaces"]
            sims[idx]["paths"] = data_reader(data)

            # Initial column placements.
            sims[idx]["cols"] = list(range(len(sims[idx]["intf"])))

            # store their intfa and intfb for assertion
            intfa.append(sims[idx]["intf"][0])
            intfb.append(sims[idx]["intf"][-1])
            intfs += sims[idx]["intf"]

    # Check that sim interfaces start and end the same place.
    assert len(set(intfa)) == 1
    assert len(set(intfb)) == 1

    # create col_order, i in (0, 1) are the 0^- and 0^+ ensembles
    tot_intfs = list(sorted(set(intfs)))
    col_order = {intf: i for i, intf in enumerate(tot_intfs[1:-1], 2)}

    # write combo toml
    with open(f"{out}.toml", "wb") as f:
        tomli_w.dump({"simulation":{"interfaces": tot_intfs}}, f)

    # merge the sim results together with the combo col
    scramble_l = []
    with open(f"{out}.txt", "w") as write:
        for idx in sims.keys():
            sim_intfs = sims[idx]["intf"]
            cols = list(range(len(sims[idx]["intf"])))
            recol = [0, 1] + [col_order[i] for i in sim_intfs[1:-1]]
            col_dic = {i: j for i, j in zip(cols, recol)}

            for line_num, path in enumerate(sims[idx]["paths"]):
                if line_num < skip:
                    continue
                frac, weig = [], []

                # create new pdic with combo cols:
                pcols = {col_dic[i]:value for i, value in path["cols"].items()}

                # fill up with weights or "----"
                for col in range(len(tot_intfs)):
                    if col in pcols:
                        frac.append(pcols[col][0])
                        weig.append(pcols[col][1])
                    else:
                        frac.append("----")
                        weig.append("----")
                string = ""
                string += f"\t{path['pn']}\t"
                string += f"{path['len']}" + "\t"
                string += f"{path['max_op']}" + "\t"
                if scramble:
                    scramble_l.append(
                        string + "\t".join(frac) + "\t" + "\t".join(weig) + "\t\n"
                    )
                else:
                    write.write(
                        string + "\t".join(frac) + "\t" + "\t".join(weig) + "\t\n"
                    )

        if scramble_l:
            np.random.seed(0)
            np.random.shuffle(scramble_l)
            for line in scramble_l:
                write.write(line)

    print("The following command can now be ran:")
    print("inft wham -toml combo.toml -data combo.txt -folder combo")
