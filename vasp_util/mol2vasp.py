# python3.X
# python mol2vasp.py ~.mol > POSCAR.vasp
import sys
import numpy as np

files = [sys.argv[1]]

for f in range(len(files)):
    o = enumerate(open(files[f], 'r'))
    NIONS = 0
    for idx, line in o:
        tmp = line.split('"')
        if 'data' in tmp:
            ref_idx = tmp.index('data')
            data = tmp[ref_idx+2].split('\\n')
            for i in range(len(data)):
                if NIONS == 0 and i > 2:
                    NIONS = int(data[i].split()[0])
                    coords_dic = {}
                    coords = []
                elif 3 < i  and i <= 3+NIONS:
                    x, y, z, ele = data[i].split()[:4]
                    if not(ele in coords_dic.keys()):
                        coords_dic[ele] = []
                    xyz = list(map(float,[x,y,z]))
                    coords_dic[ele].append(xyz)
                    coords.append(xyz)
    # Calculate the geometric center of the molecule
    # Get the longest distance between an atom and the center
    # Set the lattice vector (10 angstrom of vacuum)
    coords = np.array(coords)
    com = np.mean(coords,axis=0)
    dist = 0
    for i in range(len(coords)):
       r = np.linalg.norm(coords[i]-com)
       print(com, coords[i], r, dist, r > dist)
       if r > dist: dist = r
    lat = 2*dist+10
    trans_vec = np.array([lat/2]*3) - com

    # Generate POSCAR
    dat = 'structure from ' + files[f] + '\n'
    dat += ' 1.0\n'
    dat += f'{lat:^{20}.8f}{0:^{20}.8f}{0:^{20}.8f}\n\
{0.0:^{20}.8f}{lat:^{20}.8f}{0:^{20}.8f}\n\
{0:^{20}.8f}{0:^{20}.8f}{lat:^{20}.8f}\n'   
    types = list(coords_dic.keys())
    for ele in range(len(types)):
        dat += '  ' + types[ele]
    dat += '\n'
    for ele in range(len(types)):
        dat += '  ' + str(len(coords_dic[types[ele]]))
    dat += '\nCartesian\n'
    for ele in range(len(types)):
        for x in range(len(coords_dic[types[ele]])):
           dat += '   '.join(list(map(str,np.array(coords_dic[types[ele]][x]) + trans_vec))) + '\n'

    #pos = open(sys.argv[2]+'/POSCAR_'+files[f].split('/')[-1].split('.mol')[0],'w')
    #pos.write(dat)
    #pos.close()
    print(dat)
