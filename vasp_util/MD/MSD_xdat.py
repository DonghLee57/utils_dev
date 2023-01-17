import sys
import numpy as np

def main():
    atoms = []
    step  = 10
    n = 0

    lat = np.zeros((3,3))
    o = enumerate(open(sys.argv[1],'r'))
    next(o);next(o);
    for i in range(3):
        lat[i] = list(map(float,next(o)[1].split()))
    next(o);
    ntypes = list(map(int,next(o)[1].split()))
    NIONS = sum(ntypes)
    if len(atoms) == 0: atoms = range(NIONS)
    atom_MSD = np.array([[0.0]*len(atoms)])
    ref_coords = np.zeros((len(atoms),3))
    old_coords = np.zeros((len(atoms),3))
    new_coords = np.zeros((len(atoms),3))
    displacement = np.zeros((len(atoms),3))

    # Reference frame 0
    next(o);
    for i in range(NIONS):
        line = next(o)
        if i in atoms:
            ref_coords[atoms.index(i)] = np.matmul(list(map(float,line[1].split())),lat)
    old_coords = np.copy(ref_coords)

    # frame every "step"
    fr = 1
    for idx, line in o:
        if line[0]=='D':
            idx,line = next(o)
            fr +=1
            n = 0
        if fr % step == 0:
            if n in atoms:
                tmp = np.matmul(list(map(float,line.split())),lat)
                res = get_DIST(lat, old_coords[atoms.index(n)], tmp)
                displacement[atoms.index(n)] += res[1]
                new_coords[atoms.index(n)] = old_coords[atoms.index(n)] + res[1]
            n += 1
            if n == NIONS:
                old_coords = np.copy(new_coords)
                atom_MSD = np.linalg.norm(displacement,axis=1)**2
                p_line = f'{fr:<16d}'
                for i in range(len(atom_MSD)):
                    p_line += f'{atom_MSD[i]:^20.12f}'
                print(p_line)

def get_DIST(lat, pos_i, pos_j):
    #pos_X : Cartesian coordinates
    dist = 999
    for a in range(-1,2):
        for b in range(-1,2):
            for c in range(-1,2):
                tmp = pos_j + a*lat[0] + b*lat[1] + c*lat[2]
                eval_dist = np.linalg.norm(pos_i-tmp)
                if eval_dist < dist:
                    dist = eval_dist
                    vec = tmp-pos_i
    return dist, vec

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
