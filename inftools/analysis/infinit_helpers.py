import numpy as np
import pathlib as pl

def set_default_infinit(config):
    interfaces = config["simulation"]["interfaces"]
    assert len(interfaces) > 1, "Define some interfaces!"
    config["simulation"]["interfaces"] = [interfaces[0], interfaces[-1]]

    steps_per_iter = config["infinit"]["steps_per_iter"]
    config["infinit"]["steps_per_iter"] = steps_per_iter

    cstep = config["infinit"].get("cstep", 0)
    config["infinit"]["cstep"] = cstep
    assert cstep >= 0, "Cant restart from negative cstep??"
    assert cstep <= len(steps_per_iter[cstep:]), "Nothing to do"
    if cstep > 0:
        print(f"Restarting infinit from iteration {cstep}.")

    return config["infinit"]

def run_infretis(config, steps):
    config["simulation"]["steps"] = steps
    scheduler(config)


def get_pcross_wham(ifile, lambda_interfaces, lamres, nskip):
    # the Cxy values of [0+] are stored in the i0plus-th
    # column (first coulumn is counted as column nr 0)
    i0plus = 4

    lambdaA = lambda_interfaces[0]
    lambdaB = lambda_interfaces[-1]
    # number of interfaces
    nintf = len(lambda_interfaces)
    # number of plus-ensembles [i+]
    nplus_ens = nintf - 1
    # The sum of fractional sampling occurrences for each [i+] ensemble
    eta = np.zeros(nplus_ens)
    # Generate a list of values
    lambda_values = np.array([
        i * lamres
        for i in range(round(lambdaA / lamres), round(lambdaB / lamres) + 1)
    ])
    # v_alpha: to become the total crossing probability based on WHAM
    v_alpha = np.zeros(len(lambda_values))
    # This entry is not changed anymore
    v_alpha[0] = 1


    # Open the file
    with open(ifile) as file:
        # Initialize the matrix
        matrix = []
        # Loop through each line in the file
        for line in file:
            # Ignore comment lines
            if line.startswith("#"):
                continue
            # Remove any trailing whitespace or newlines
            line = line.strip()
            # Split the line into a list of values
            values = line.split()
            # Replace any "----" values with 0.0
            values = [float(x) if x != "----" else 0.0 for x in values]
            # Add the row to the matrix
            matrix.append(values)

    # delete the first nskip entries
    del matrix[:nskip]

    # Format of data-file is not yet unweighted with HA-weights
    # column for [0-] ensemble
    i0min = i0plus - 1
    # sum of Pxy for all y being [0-], [0+], [1+] etc
    sumPxy = [0.0] * nintf
    # sum of Pxy after weighting with inverse HA-weights
    sumPxy_afterw = [0.0] * nintf
    for x in matrix:
        for y in range(nintf):
            # index of value that requires unweighting
            y1 = i0min + y
            # index of HA-weight
            y2 = y1 + nintf
            if x[y2] > 0:  # non-zero weight
                sumPxy[y] += x[y1]  # sum before unweighting
                x[y1] /= x[y2]
                sumPxy_afterw[y] += x[y1]  # sum after unweighting
            elif x[y1] > 0:
                print("Division by zero for path x=", x)
                exit()
            # sumPxy and suminvw remain the same.
            # if x[y1]==x[y2]==0, giving 0/0, x[y1] remains zero.
            # Next:
            # we used to store the HA-weights
            # store the running sum of sumPxy in the matrix where
            # running SumPxy values are needed for
            # computing the running WHAM averages
            # HA-weights are no longer needed after this step, while
            x[y2] = sumPxy[y]
    # loop over x completed
    # dvide each column by the average of the inverse HA-weight
    # This gives more comparible eta[y] values between ensembles.
    # For instance if [0+] is based on shooting and [i+] is based on WF
    # This allows to use the standard WHAM procedure
    for y in range(nintf):
        AvinvwHA = sumPxy_afterw[y] / sumPxy[y]
        y1 = i0min + y  # index of [0-], [0+], [1+] etc
        for x in matrix:
            x[y1] /= AvinvwHA

    # This function gives back the WHAM based crossing
    # probabilities at the interfaces and the Q-factor
    # to normailze the v_alpha vector
    def WHAM_PQ(npe, interf, res, eta, v_alpha):
        P = [0.0] * npe  # Crossin probability at lambda_0, lambda_1,
        #  .... lambda_{n-1}: P_A(\lambda_i|\lambda_0) based on WHAM
        Q = [0.0] * npe  # The Q-factor for Whamming (JCTC 2015 Lervik et al)
        invQ = [0.0] * npe  # The inverse Q-factor
        # Initial values
        P[0] = 1.0
        invQ[0] = eta[0]
        if invQ[0] == 0:
            return P, Q
        Q[0] = 1 / invQ[0]
        # Solve other values using recursive relations
        lambdaA = interf[0]
        for i in range(1, npe):  # Loop over all nterfaces after [0+]
            lambda_i = interf[i]
            alpha = round(
                (lambda_i - lambdaA) / res
            )  # the alpha index corresponding to lambda_i
            P[i] = v_alpha[alpha] * Q[i - 1]
            if P[i] == 0:
                return P, Q
            invQ[i] = invQ[i - 1] + (eta[i] / P[i])
            Q[i] = 1 / invQ[i]
        # add the final value of lambda_B
        P.append(v_alpha[-1] * Q[nplus_ens - 1])
        return P, Q

    # Loop over all paths x in the matrix to:
    # Compute eta[i] for all interfaces except the last lambda_B
    # Create the v_alpha() vector. Yet, without the proper normalization
    # Create the local crossing probabilities p_loc[.... ]
    # which is a matrix showing the local crossing probabilities
    # for [0+], [1+] etc
    # Create running average of local crossing probability
    # P_A(\lambda_{i+1}|\lambda_i): run_av_ploc
    # The structure of the data file is the following
    # x[0]=path-index, x[1]=path length, x[2]=lambda-max
    for x in matrix:
        lambdamax = x[2]
        for i in range(nplus_ens):
            Cxy_index = i0plus + i
            Cxy = x[Cxy_index]
            eta[i] += Cxy  # increase eta[i]
            # Determine lower and upper bound for
            # increasing v_alpha- and p_loc values
            lambda_i = lambda_interfaces[i]
            alpha_max = int(np.floor((lambdamax - lambdaA) / lamres))
            alpha_min = round((lambda_i - lambdaA) / lamres)
            # Note: lambda_i-lambdaA)/lamres is an integer as
            # lambda_i and lambdaA should be commensurate with lamres
            if alpha_max > len(v_alpha) - 1:
                alpha_max = len(v_alpha) - 1  # -1 as we start counting from 0
            alpha_min += 1  # v(alpha) at the interface lambda_1, lambda_2
            # etc are determined by the previous [0+], [1+] etc
            for alpha in range(alpha_min, alpha_max + 1):
                v_alpha[alpha] += Cxy

    # Normalize the p_loc arrays such that each local
    # crossing probability of ensemble [i+] starts with 1 at lambda_i
    _, Q = WHAM_PQ(nplus_ens, lambda_interfaces, lamres, eta, v_alpha)

    # we need now a loop ovber all alpha values and determine K(alpha)
    # It is however faster to loop over K(alpha) indexes and split the
    # alpha-loop into parts that have the same K(alpha) value
    for Kalpha in range(nplus_ens):  # Loop over all pluse-ensembles after [0+]
        lambda_i = lambda_interfaces[Kalpha]
        lambda_next = lambda_interfaces[Kalpha + 1]
        alpha_max = round((lambda_next - lambdaA) / lamres)
        # Note again (lambda_next-lambdaA)/lamres is an integer
        # as interfaces should be commensurate with lamres resolution
        alpha_min = (
            round((lambda_i - lambdaA) / lamres) + 1
        )  # +1 because K(lambda) refers to index K with lambda_K < lambda
        # (not less or equal)
        for alpha in range(alpha_min, alpha_max + 1):
            v_alpha[alpha] *= Q[Kalpha]
    # v_alpha now represents the total crosing probability based on WHAM
    return lambda_values, v_alpha


def print_logo(step: int = 0):
    from rich.console import Console
    from rich.text import Text

    console = Console()
    art = Text()
    # add some initial states, a [0+] path and [0-] path,
    # some A-B paths, etcetera etcetera
    if step >= len(PHRASES) or step == -1:
        step = np.random.choice(len(PHRASES))
    str1,str2 = PHRASES[step]
    #art.append("\n")
    art.append("         o             ∞       ____          \n", style="bright_blue")
    art.append("__________\________o__________/____\_________\n", style="bright_blue")
    art.append("           \      /          o      o        \n", style="bright_blue")
    art.append("            \____/                           \n", style="bright_blue")
    art.append(str1,style="bright_cyan")
    art.append(f"{'_'*(45-len(str1))}\n", style="bright_blue")
    art.append("_____________________________________________\n", style="bright_blue")
    art.append(" _          __  _         _  _               \n", style="bold light_cyan")
    art.append("(_) _ __   / _|(_) _ __  (_)| |_             \n", style="bold bright_yellow")
    art.append("| || '_ \ | |_ | || '_ \ | || __|            \n", style="bold bright_magenta")
    art.append("| || | | ||  _|| || | | || || |_             \n", style="bold bright_green")
    art.append("|_||_| |_||_|  |_||_| |_||_| \__|            \n", style="bold white")
    art.append("______________\______________________________\n", style="bright_blue")
    art.append("   ∞           o                o            \n", style="bright_blue")
    #art.append("             o                    ",style="bright_blue")
    art.append(f"{str2:>45}", style="bright_cyan")
    console.print(art)

class LightLogger:
    def __init__(self, fname):
        self.fname = str(pl.Path(fname).resolve())

    def log(self, msg):
        with open(self.fname, "a") as wfile:
            wfile.write(msg + "\n")
