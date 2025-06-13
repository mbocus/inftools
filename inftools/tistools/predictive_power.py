import typer
import tqdm
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Annotated, Optional

def predictive_power(
    lambda_c: Annotated[float, typer.Option("-lambda_c", help="The op value of the interface at which the CV is evaluated")],
    lambda_r: Annotated[float, typer.Option("-lambda_r", help="The op value of the reactive interface")],
    fn_weights: Annotated[str, typer.Option("-fn_weights", help="The path to the paths weight list")],
    cv_col: Annotated[int, typer.Option("-cv_col", help="The column of the order.txt file at which the cv can be found")] = None,
    op_col: Annotated[int, typer.Option("-op_col", help="The column of the order.txt file at which the infretis op can be found")] = 1,
    cv_file: Annotated[str, typer.Option("-cv_file", help="The path to a Python file that defines a 'calc_cv(atoms)' function, where atoms is an ASE Atoms object")] = None,
    cv_bins: Annotated[int, typer.Option("-cv_bins", help="The number of bins to build the cv histograms")] = 10,
    savgol: Annotated[bool, typer.Option("-savgol/-no-savgol", help="Use a Savitzky-Golay filter for binless analysis. If true, cv_bins is only used for plotting", show_default=True)] = True,
    path_load: Annotated[str, typer.Option("-path_load", help="The path to the infretis load folder")] = "load",
):
    """
    Compute the predictive power of a certain collective variable (cv) given a first interface placed at lambda_c
    (crossed interface) and a second interface placed at lambda_r (reactive interface). The weights file can 
    be obtained with the command:
        inft get_path_weights -data infretis_data.txt -toml infretis.toml -nskip 50 -out path_weights.txt
    The new cv can be either precomputed as column of the order.txt file or computed on the fly with a custom function.
    The savgol option uses a Savitzky-Golay filter to reduce the dependency on the bins, see https://doi.org/10.1021/acs.jctc.5c00054. It relies on default options,
    which can be customized within this script.

    Current limitations:
    - can only work on a single CV
    - the cv calculation is ASE centric, it expects to work with Atoms
    - the function only accepts one value for lambda c and one for lambda r
    """

    assert (cv_col is None) != (cv_file is None), "Please provide either cv_col or cv_file, but not both"
    if cv_col is not None and cv_col == op_col:
        raise ValueError("cv_col and op_col must be different. The cv should be different from the original order parameter.")

    assert lambda_c < lambda_r, "lambda_c should be smaller than lambda_r!"

    # initialize weights and cv values for reactive and unreactive paths
    cvu, wu = [], []
    cvr, wr = [], []

    weights = np.loadtxt(fn_weights) # the weights file, also contains the max_op and path ID
    for line in tqdm.tqdm(weights):
        path_id = int(line[0])
        path_dir = Path(path_load) / f"{path_id}"
        max_op = line[1]
        weight = line[2]

        # load the order.txt file of the path
        order = np.loadtxt(path_dir / "order.txt")[1:]
        for cvs in order:
            if cvs[op_col] > lambda_c: # check where lambda_c is crossed
                if cv_col is not None:
                    cv = cvs[cv_col]

                else: # compute the cv on the fly
                    import importlib.util
                    from ase.io import read
                    path_cv_file = Path(cv_file)
                    # make sure the python file exists
                    if not path_cv_file.is_file() or path_cv_file.suffix != ".py":
                        typer.echo(f"Error: {cv_file} does not exist or it is not a .py file!")
                        raise typer.Exit(code=1)
                    
                    # import the cv module
                    module_name = path_cv_file.name
                    spec = importlib.util.spec_from_file_location(module_name, cv_file)
                    if spec is None:
                        typer.echo(f"Error: Could not load module from {cv_file}.")
                        raise typer.Exit(code=1)
                    module = importlib.util.module_from_spec(spec)
                    sys.modules[module_name] = module
                    try:
                        spec.loader.exec_module(module)
                    except Exception as e:
                        typer.echo(f"Error importing module: {e}")
                        raise typer.Exit(code=1)
                    
                    # make sure that calc_cv exists and its callable
                    if not hasattr(module, 'calc_cv') or not callable(module.calc_cv):
                        typer.echo(f"Error: Module '{module_name}' does not contain a callable 'calc_cv' function.")
                        raise typer.Exit(code=1)
                    
                    # identify the structure corresponding to the order value
                    frame_index = cvs[0]
                    with open(path_dir / "traj.txt", 'r') as traj:
                        lines = (line for line in traj if not line.lstrip().startswith('#'))
                        for i, line in enumerate(lines):
                            if i == frame_index:
                                fn_traj = line.split()[1]
                                traj_index = int(line.split()[2])
                                break
                    # load the atoms and compute the cv
                    atoms = read(path_dir / f"accepted/{fn_traj}", index=traj_index)
                    cv = module.calc_cv(atoms)
                
                # check if the path is reactive or not
                if max_op > lambda_r:
                    cvr.append(cv)
                    wr.append(weight)
                else:
                    cvu.append(cv)
                    wu.append(weight)
                break # no need to continue
    
    cvr, wr = np.array(cvr), np.array(wr)
    cvu, wu = np.array(cvu), np.array(wu)

    W = np.sum(wr) + np.sum(wu)
    wr /= W
    wu /= W
    Pr = np.sum(wr)
    print("Crossing probability at lambda_r:", Pr)

    if savgol:
        from scipy.signal import savgol_filter
        import pandas as pd

        ### Filter settings ###
        window_divisor = 16 # The filter will consider 1 / window_divisor of the data at a time
        pol_deg = 2 # degree of the polynomial interpolation
        ext_factor = 1.5 # increases the cv range by this factor to account for the plateaus at the beginning and end
        n_data_grid = 2000 # the number of grid points in the original cv range
        #######################

        n_grid = int(n_data_grid * ext_factor)
        n_start = int((n_grid - n_data_grid) / 2)
        n_end = int(n_start + n_data_grid + 1 )

        df_r = pd.DataFrame({"cv": cvr, "w": wr})
        df_u = pd.DataFrame({"cv": cvu, "w": wu})
        # reorder the cv and weights such that the cv values are in ascending order
        df_r.sort_values("cv", ascending=True, inplace=True)
        df_u.sort_values("cv", ascending=True, inplace=True)
        # compute the cumulative sum of the weights
        df_r["w"] = df_r["w"].cumsum()
        df_u["w"] = df_u["w"].cumsum()
        # drop the cv duplicates (if any)
        df_r.drop_duplicates(subset="cv", keep="last", inplace=True)
        df_u.drop_duplicates(subset="cv", keep="last", inplace=True)

        # prepare extended cv grid for savgol interpolation
        cvmin, cvmax = min([df_r["cv"].iloc[0], df_u["cv"].iloc[0]]), max([df_r["cv"].iloc[-1], df_u["cv"].iloc[-1]])
        cvrange = cvmax - cvmin
        cvmid = cvmin + 0.5 * cvrange
        ext_cvrange = cvrange * ext_factor
        ext_cvmin, ext_cvmax = cvmid - 0.5 * ext_cvrange, cvmid + 0.5 * ext_cvrange
        ext_dcv = ext_cvrange / n_grid
        ext_cv = np.arange(ext_cvmin, ext_cvmax, ext_dcv)

        # linear interpolation between datapoints
        wr_linear = np.zeros(len(ext_cv))
        wu_linear = np.zeros(len(ext_cv))
        wr_linear[n_start:n_end] += np.interp(ext_cv[n_start:n_end], df_r["cv"], df_r["w"])
        wu_linear[n_start:n_end] += np.interp(ext_cv[n_start:n_end], df_u["cv"], df_u["w"])
        if n_grid > n_data_grid:
            wr_linear[n_end:] += wr_linear[n_end-1]
            wu_linear[n_end:] += wu_linear[n_end-1]

        # Savitzky-Golay interpolation (note that we ask for the first derivative of the curve)
        window_length = int(n_data_grid / window_divisor)
        if (window_length % 2) == 0:
            window_length += 1
        wr_savgol = savgol_filter(wr_linear, window_length, pol_deg, deriv=1, delta=ext_dcv)
        wu_savgol = savgol_filter(wu_linear, window_length, pol_deg, deriv=1, delta=ext_dcv)

        # round such that the smallest value is two orders of magnitude smaller than Pr
        decimals = int(abs(np.floor(np.log10(abs(Pr)))) + 2)
        wr_savgol = np.round(wr_savgol, decimals)
        wu_savgol = np.round(wu_savgol, decimals)

        wr_savgol[wr_savgol==0] = np.nan
        wu_savgol[wu_savgol==0] = np.nan

        # predictive power
        overlap = (1 / Pr) * wu_savgol * wr_savgol * ext_dcv / (wu_savgol + wr_savgol)
        pp = 1 - np.nansum(overlap)
        print("Computed predictive power:", pp)

        # plot results
        fig, axs = plt.subplots(1, 2, figsize=(8,4))
        bins = np.linspace(ext_cv[0], ext_cv[-1], cv_bins)
        # reactive paths
        axs[0].plot(ext_cv, wr_savgol, label="SG interpolation")
        counts, bins, patches = axs[0].hist(cvr, bins=bins, weights=wr, alpha=0.5)
        # normalize histogram height based on the bin width
        bin_widths = np.diff(bins)
        for count, patch, width in zip(counts, patches, bin_widths):
            patch.set_height(count / width)
        axs[0].set_title("Reactive")
        # unreactive paths
        axs[1].plot(ext_cv, wu_savgol, label="SG interpolation")
        counts, bins, patches = axs[1].hist(cvu, bins=bins, weights=wu, alpha=0.5)
        # normalize histogram height based on the bin width
        bin_widths = np.diff(bins)
        for count, patch, width in zip(counts, patches, bin_widths):
            patch.set_height(count / width)
        axs[1].set_title("Unreactive")
        for ax in axs.flatten():
            ax.legend()
            ax.set_xlabel("CV")
        plt.show()

    else:
        W = np.sum(wr) + np.sum(wu)

        # simply compute histograms
        cvtot = np.concatenate((cvr,cvu))
        cvmin, cvmax = np.min(cvtot), np.max(cvtot)
        bins = np.linspace(cvmin, cvmax, cv_bins)
        hr, hr_bin_edges = np.histogram(cvr, bins=bins, weights=wr)
        hu, hu_bin_edges = np.histogram(cvu, bins=bins, weights=wu)

        # predictive power
        overlap = (1 / Pr) * hr * hu / (hr + hu)
        pp = 1 - np.nansum(overlap)
        print("Computed predictive power:", pp)

        # plot results
        fig, ax = plt.subplots()
        ax.stairs(hr, hr_bin_edges, label="reactive")
        ax.stairs(hu, hu_bin_edges, label="unreactive")
        ax.set_yscale('log')
        ax.set_xlabel("CV")
        ax.legend()
        plt.show()
