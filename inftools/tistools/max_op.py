import argparse

import numpy as np
import matplotlib.pyplot as plt
import tomli

COLS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
        '#bcbd22', '#17becf']

def infretis_data_reader(input):
    with open(input, 'r') as read:
        _, head, _ = (read.readline(), read.readline(), read.readline())
        num_ens = len(head.rstrip().split()) - 5
        ens_dic = {i : {'op': [], 'haw': [], 'w': []} for i in range(num_ens)}
        for line in read:
            rip = line.rstrip().split()[3:]
            pinfo = line.rstrip().split()[:3]
            for idx, (w, haw) in enumerate(zip(rip[:num_ens], rip[num_ens:])):
                if '----' in (w, haw):
                    continue
                ens_dic[idx]['op'].append(float(pinfo[2]))
                ens_dic[idx]['w'].append(float(w))
                ens_dic[idx]['haw'].append(float(haw))
    return ens_dic

def plot_max_op(arguments):
    parser = argparse.ArgumentParser(
        description="Generate initial paths for an infretis \
                    simulation using paths from an earlier infretis \
                    simulation"
    )

    parser.add_argument(
        "-data",
        help="The path to the infretis_data.txt data file",
    )
    parser.add_argument(
        "-toml",
        help="The .toml input file for reading the interfaces\
                (e.g. ../iretis0/infretis.toml)",
    )
    parser.add_argument(
        "-ensp",
        help="the plus ensembles",
        type=int,
    )

    args = parser.parse_args(arguments)
    # read interfaces from .toml file
    with open(args.toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    intfs = toml_dict["simulation"]["interfaces"]

    plt.axhline(intfs[0], color=f'k')
    plt.axhline(intfs[-1], color=f'k')
    plt.axhline(intfs[args.ensp], color=f'C0')
    plt.axhline(intfs[args.ensp+1], color=f'C1')
    print('# interface', intfs[args.ensp])

    ens_dic = infretis_data_reader(args.data)
    data = ens_dic[args.ensp+1]
    # print('len', len(data['op']))
    # we = np.array(data['w'])/np.array(data['haw'])
    # we = we/max(we)
    for i in range(2):
        we = np.array(data['w'])/np.array(data['haw'])
        we = we/max(we)
        maxidx = np.argmax(we)
        data['op'].pop(maxidx)
        data['haw'].pop(maxidx)
        data['w'].pop(maxidx)
    we = np.array(data['w'])/np.array(data['haw'])
    we = (we/max(we))**2
    # print(sorted(we))
    for i, j in enumerate(data['op']):
        print(i, j, data['w'][i], data['haw'][i])

    plt.scatter(list(range(len(data['op']))), data['op'], color=f'C0')#, alpha=we)

    # for idx, (key, item) in enumerate(ens_dic.items()):
    #     if idx == 0 or idx - 1 not in args.ensp:
    #         print(idx - 1, args.ensp)
    #         continue
    #     # if idx - 1 not in args.ensp:
    #     plt.scatter(list(range(len(item['op']))), item['op'], color=f'C{idx-1%8}')
    #     # break
    plt.show()
        # print(key, len(item['op']))

    # print(ens_dic)
    # with open(args.data, 'r') as read:
    #     _, head, _ = (read.readline(), read.readline(), read.readline())
    #     num_ens = len(head.rstrip().split()) - 3
    #     ens_dic = {i : {'op': [], 'haw': []} for i in range(num_ens)}
    #     for line in read:
    #         rip = line.rstrip().split()[3:]
    #         for idx, (op, haw) in enumerate(zip(rip[:num_ens], rip[num_ens:])):
    #             if '----' in (op, haw):
    #                 continue
    #             print('whaid', (op, haw), '-' in (op, haw))
    #             ens_dic[idx]['op'].append(float(op))
    #             ens_dic[idx]['haw'].append(float(haw))








