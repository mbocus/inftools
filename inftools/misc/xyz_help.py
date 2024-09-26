def calc_center(center, xyzs, cell):
    xyzs_p = []
    for xyz in xyzs:
        xyz_p = [xyz[0]]
        for dim, center_i in zip(xyz[1:], center):
            dim_c = dim - center_i
            if dim_c < -0.5 * cell:
                mult = 1
            elif dim_c > 0.5 * cell:
                mult = -1
            else:
                mult = 0
            xyz_p.append(f"{dim_c+cell*mult:10.8f}")
        xyzs_p.append("\t".join(xyz_p) + "\n")
    return xyzs_p
