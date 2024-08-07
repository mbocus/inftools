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

