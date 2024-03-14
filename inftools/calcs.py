from os import walk
import matplotlib.pyplot as plt
import numpy as np

from inftools.plotting.max_op import infretis_data_reader


#def read_infdata(inp):
#    dic = {'pn': [], 'len': [], 'max_op': []}
#    ens_dic = {}
#
#    with open(inp, 'r') as read:
#        for line in read:
#            if '#' in line:
#                continue
#            rip = line.rstrip().split()
#            dic['pn'].append(int(rip[0]))
#            dic['len'].append(int(rip[1]))
#            dic['max_op'].append(float(rip[2]))
#            ensnu = len(rip[3:][::2])
#            ensv, ensw = rip[3:][:ensnu], rip[3:][ensnu:]
#            for idx, (weight, haweight) in enumerate(zip(ensv, ensw)):
#                if '--' in weight + haweight:
#                    continue
#                else:
#                    ens_n = f'{idx:03}'
#                    if ens_n not in ens_dic:
#                        ens_dic[ens_n] = []
#                    ens_dic[ens_n].append([float(weight), float(haweight), float(rip[2])])
#
#    # print(sorted(ens_dic.keys()))
#    # print(ens_dic['000'])
#    # exit('tiger')
#    return dic, ens_dic



def calc_rate(arguments):
    import argparse

    parser = argparse.ArgumentParser(
        description="Estimate interfaces from Pcross.txt"
    )
    parser.add_argument("-wham", help="The wham folder")
    parser.add_argument("-dt", type=float, help="timestep")
    parser.add_argument("-sc", type=float, help="subcycle")
    parser.add_argument("-unit", help="defaults: cp2k: 'fs', gromacs: 'ps'")
    parser.add_argument("-plot", help="specific value(s) vs runavg, in s-1.", nargs="+",)

    args = parser.parse_args(arguments)

    run_rate = np.loadtxt(args.wham + '/runav_rate.txt')
    run_flux = np.loadtxt(args.wham + '/runav_flux.txt')
    run_pcro = np.loadtxt(args.wham + '/runav_Pcross.txt')

    run_rate[:, 3] = run_rate[:, 3]/(args.dt*args.sc)
    run_flux[:, 2] = run_flux[:, 2]/(args.dt*args.sc)
    print('calculated rate:', f'{run_rate[-1, 3]:.3e}', f"{args.unit}^-1")
    if args.unit == 'fs':
        run_rate[:, 3] = run_rate[:, 3]*10**(15)
        run_flux[:, 2] = run_flux[:, 2]*10**(15)
        print('calculated rate:', f'{run_rate[-1, 3]:.3e}', "s^-1")
        print('calculated flux:', f'{run_flux[-1, 2]:.3e}', "s^-1")
        print('calculated p_ab:', f'{run_pcro[-1, 2]:.3e}', "s^-1")
        print('calculated rate:', f'{run_flux[-1, 2]*run_pcro[-1, 2]:.3e}', "s^-1")

        print(run_rate[-1, 3], "s^-1")

    if args.plot:
        import scienceplots
        plt.style.use(["science"])
        plt.plot(run_rate[:, 0], run_rate[:, 3], label=r'$\infty$RETIS',
                 color='#1f77b4', linewidth=2.0)
        for i in args.plot:
            minval = min(run_rate[-1, 3], float(i))
            maxval = max(run_rate[-1, 3], float(i))
            print(i, "s^-1", ', difference, a factor of:', maxval/minval)
            plt.axhline(float(i), color='k', label='Experimental', linewidth=2.0)
        plt.legend(frameon=False)
        plt.yscale('log')
        plt.xlabel('Number of accepted paths')
        plt.ylabel(r'Rate constant [s$^{-1}$]')
        plt.savefig('co2_rate.pdf')
        # plt.show()



def calc_ensdata(arguments):
    import argparse
    import tomli

    parser = argparse.ArgumentParser(
        description="Estimate interfaces from Pcross.txt"
    )
    parser.add_argument("-data", help="The infretis_data.txt file", default='infretis_data.txt')
    parser.add_argument("-plot", action="store_true", help="Plot")
    parser.add_argument(
        "-toml", help="The .toml input file defining the orderparameter", default='infretis.toml'
    )

    args = parser.parse_args(arguments)
    with open(args.toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    interfaces = toml_dict["simulation"]["interfaces"]

    # dic, ens_dic = read_infdata(args.data)
    ens_dic2 = infretis_data_reader(args.data)
    # keysa = list(sorted(ens_dic.keys()))
    # keysb = list(sorted(ens_dic2.keys()))
    # print('gori a', ens_dic.keys())
    # print('gori b', ens_dic2.keys())
    # for keya, keyb in zip(keysa, keysb):
    #     print(keya, keyb)
    #     print(len(ens_dic[keya]))
    #     print(len(ens_dic2[keyb]['haw']))
    # exit('uhuh')
    enss = sorted(ens_dic2.keys())
    num = [len(i['op']) for i in ens_dic2.values()]
    # for intf, num0 in zip(interfaces[:-1], num[1:]):
    #     print(intf, num0)

    # for idx, ens in enumerate(ens_dic.values()):
    pcross = 1.
    pcrosses = []
    for ens in enss:
        # skip zero minus ensemble
        data = ens_dic2[ens]
        if ens == 0:
            # print(data)
            # exit('uhuh')
            continue
        avg_up = 0
        avg_dw = 0
        for op, haw, w in zip(data['op'], data['haw'], data['w']):
            if op > interfaces[ens]:
                avg_up += 1*w/haw
            avg_dw += w/haw

        # plt.scatter(list(range(len([i[2] for i in data]))), [i[2] for i in data])
        # plt.axhline(interfaces[idx-1], color='C0')
        # plt.axhline(interfaces[idx], color='C1')
        # plt.show()
        # print([i[2] for i in data])
        # print(pn[-1], pn[2] > interfaces[idx+1], interfaces[idx+1])
        # exit('noouh')
        print(ens, avg_up, avg_dw)
        print('whada', avg_up/avg_dw)
        pcross *= avg_up/avg_dw
        pcrosses.append(avg_up/avg_dw)
    print('blud', pcross, len(pcrosses))
    # exit('uz')
        # exit('noouh')


        # print('a', ens)

    print('whad', pcrosses)
    if args.plot:
        plt.scatter(interfaces[:-1], num[1:],)
        plt.plot(interfaces[:-1], num[1:],)
        plt.show()
    for intf, num0, pc in zip(interfaces[:-1], num[1:], pcrosses):
        print(intf, num0, pc)

def get_time(arguments):
    import argparse
    from datetime import datetime

    parser = argparse.ArgumentParser(
        description="Estimate interfaces from Pcross.txt"
    )
    parser.add_argument("-i", help="The log file", default='sim.log')
    parser.add_argument("-data", help="The infretis_data.txt file", default='infretis_data.txt')
    parser.add_argument("-nskip", type=int, help="time from nskip", default=0)
    parser.add_argument("-plot", help="plot the data", action='store_true')

    args = parser.parse_args(arguments)
    pnum = 0
    if args.nskip!= 0:
        with open(args.data, 'r') as read:
            _, _, _ = read.readline(), read.readline(), read.readline()
            for idx, line in enumerate(read):
                pnum = line.rstrip().split()[0]
                if idx == args.nskip:
                    break

    starttime = 0
    times = []
    delta = 0
    skip_idx = -1
    with open(args.i, 'r') as read:
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
    tottime_d = int(np.floor(tottime))
    tottime_h = (tottime%1)*24
    t_label = f'{tottime:.02f} D / {tottime_d} D {tottime_h:.1f} H'
    print('total time:', t_label)
    skiptime = times[-1] - times[skip_idx]
    skiptime_d = int(np.floor(skiptime))
    skiptime_h = (skiptime%1)*24
    s_label = f'{skiptime:.02f} D / {skiptime_d} D {skiptime_h:.1f} H'
    if args.nskip > 0:
        print('skip time:', s_label)
    
    if args.plot:
        plt.plot(list(range(len(times))), times, label=t_label)
        if args.nskip > 0:
            plt.axvline(skip_idx, color='k', label=s_label)
        plt.legend(frameon=False)
        plt.ylabel('simulation time [days]')
        plt.show()
