def check_data(arguments):
    import numpy as np
    import matplotlib.pyplot as plt
    import tomli
    import argparse

    parser = argparse.ArgumentParser(description = """Check that all ensembles [i] contain paths
    that cross the interface of the next ensemble [i+1]. Can be used to plot the maximum OP value of all paths within each ensemble.""")
    parser.add_argument("-data", default = "infretis_data.txt", help = "the infretis_data.txt file")
    parser.add_argument("-toml", default = "infretis.toml", help = "the .toml file, used to read interfaces")
    parser.add_argument("--plot", default = False, action=argparse.BooleanOptionalAction, help = "plot the max OP for all paths in all ensembles")

    args = parser.parse_args(arguments)

    with open(args.toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)

    interfaces = toml_dict["simulation"]["interfaces"][:]
    x = np.loadtxt(args.data, dtype = str)

    N_ens = len(interfaces)

    # array where element [i,j] is the maximum OP
    # value of line i in ensemble j.
    OP_matrix = np.zeros((x.shape[0],N_ens))
    for i, xi in enumerate(x):
        x_tmp = x[i, 3:N_ens + 3]
        idx = np.where(x_tmp!="----")
        # drop load paths which are all '----'
        if len(idx[0])==0:
            continue

        OP_matrix[i, idx[0]] = xi[2].astype(float)

    # histogram bin width stuff

    if args.plot:
        f,a = plt.subplots(1, N_ens, sharey = True)

    for ens in range(N_ens):

        # Left and Middle interfaces
        ens_intf = ens
        if ens > 0:
            ens_intf -=1
        L = interfaces[ens_intf]
        M = interfaces[ens_intf + 1]

        # histogram
        path_nr = OP_matrix[:,ens]!=0
        op = OP_matrix[:,ens][path_nr]

        if args.plot:
            a[ens].axvline(L, c = 'k', lw=3)
            a[ens].axvline(M, c='k', lw=3)
            a[ens].scatter(op,x[path_nr, 0].astype(int))
            a[ens].set(xlim=(L, max(max(op), M)), xlabel = 'maximum OP', ylabel = 'path nr. in infretis_data.txt', title = f"ensemble: {ens}, intf={M}")
        if ens > 0 and max(op) < M:
            print(f"Ensemble {ens} is missing paths that cross interface {ens_intf+1} with value {M}")

    if args.plot:
        plt.show()
