from typing import Annotated
import typer

import matplotlib.pyplot as plt
import numpy as np

# from inftools.plotting.max_op import infretis_data_reader

def calc_rate(
    dt: Annotated[float, typer.Option("-dt", help="Timestep")],
    sc: Annotated[int, typer.Option("-sc", help="Subcycle")],
    unit: Annotated[str, typer.Option("-unit", help="cp2k: 'fs', gromacs: 'ps'")],
    plot: Annotated[float, typer.Option("-plot", help="specific value(s) vs runavg, in s-1.")] = None,
    wham: Annotated[str, typer.Option("-wham", help="The wham folder")] = "wham",
        ):
    """Calculates the rate with correct units. Currently unfinished."""
    run_rate = np.loadtxt(wham + '/runav_rate.txt')
    run_flux = np.loadtxt(wham + '/runav_flux.txt')
    run_pcro = np.loadtxt(wham + '/runav_Pcross.txt')

    run_rate[:, 3] = run_rate[:, 3]/(dt*sc)
    run_flux[:, 2] = run_flux[:, 2]/(dt*sc)
    print('calculated rate:', f'{run_rate[-1, 3]:.3e}', f"{unit}^-1")
    if unit == 'fs':
        run_rate[:, 3] = run_rate[:, 3]*10**(15)
        run_flux[:, 2] = run_flux[:, 2]*10**(15)
        print('calculated rate:', f'{run_rate[-1, 3]:.3e}', "s^-1")
        print('calculated flux:', f'{run_flux[-1, 2]:.3e}', "s^-1")
        print('calculated p_ab:', f'{run_pcro[-1, 2]:.3e}', "s^-1")
        print('calculated rate:', f'{run_flux[-1, 2]*run_pcro[-1, 2]:.3e}', "s^-1")

        print(run_rate[-1, 3], "s^-1")

    if plot:
        import scienceplots
        plt.style.use(["science"])
        plt.plot(run_rate[:, 0], run_rate[:, 3], label=r'$\infty$RETIS',
                 color='#1f77b4', linewidth=2.0)
        if plot is not None:
            minval = min(run_rate[-1, 3], float(i))
            maxval = max(run_rate[-1, 3], float(i))
            print(plot, "s^-1", ', difference, a factor of:', maxval/minval)
            plt.axhline(float(i), color='k', label='Experimental', linewidth=2.0)
        plt.legend(frameon=False)
        plt.yscale('log')
        plt.xlabel('Number of accepted paths')
        plt.ylabel(r'Rate constant [s$^{-1}$]')
        plt.savefig('co2_rate.pdf')

def get_time(
    log: Annotated[str, typer.Option("-log", help="The sim.log file")] = "sim.log",
    data: Annotated[str, typer.Option("-data", help="The infretis_data.txt file")] = "infretis_data.txt",
    nskip: Annotated[int, typer.Option("-nskip", help="Time from nskip")] = 0,
    plot: Annotated[bool, typer.Option("-plot", help="Time from nskip")] = False,
        ):
    """Calculates the time taken per path. Currently unfinished."""
    from datetime import datetime

    pnum = 0
    if nskip!= 0:
        with open(data, 'r') as read:
            _, _, _ = read.readline(), read.readline(), read.readline()
            for idx, line in enumerate(read):
                pnum = line.rstrip().split()[0]
                if idx == nskip:
                    break

    starttime = 0
    times = []
    delta = 0
    skip_idx = -1
    with open(log, 'r') as read:
        for idx, line in enumerate(read):
            if "submit" in line and "END" in line:
                date = line.split()[-2:]
                time = [int(i) for i in date[1].split(':')]
                date = [int(i) for i in date[0][7:].split('.')]
                starttime = datetime(date[0], date[1], date[2],
                                     time[0], time[1], time[2])
                if times:
                    delta = times[-1]

            if "date" in line:
                date = line.split()[-2:]
                time = [int(i) for i in date[1].split(':')]
                date = [int(i) for i in date[0].split('.')]
                ctime  = datetime(date[0], date[1], date[2],
                                  time[0], time[1], time[2])
                times.append((ctime - starttime).total_seconds() + delta)
            if f'p{pnum}' in line and skip_idx == -1:
                skip_idx = len(times)

    times = [i/(60*60*24) for i in times]
    tottime = times[-1] # days
    print('gorilla', tottime)
    tottime_d = int(np.floor(tottime))
    tottime_h = (tottime%1)*24
    t_label = f'{tottime:.02f} D / {tottime_d} D {tottime_h:.1f} H'
    print('total time:', t_label)
    skiptime = times[-1] - times[skip_idx]
    skiptime_d = int(np.floor(skiptime))
    skiptime_h = (skiptime%1)*24
    s_label = f'{skiptime:.02f} D / {skiptime_d} D {skiptime_h:.1f} H'
    if nskip > 0:
        print('skip time:', s_label)

    if plot:
        plt.plot(list(range(len(times))), times, label=t_label)
        if nskip > 0:
            plt.axvline(skip_idx, color='k', label=s_label)
        plt.legend(frameon=False)
        plt.xlabel('Number of paths.')
        plt.ylabel('Simulation time [days]')
        plt.show()
