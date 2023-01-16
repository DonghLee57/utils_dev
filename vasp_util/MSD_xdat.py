import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    atoms = []
    MSD = [0.0]
    step = 10

    lat = np.zeros((3,3))
    o = enumerate(open(sys.argv[1],'r'))
    next(o);next(o);
    for i in range(3):
        lat[i] = list(map(float,next(o)[1].split()))
    next(o);
    ntypes = list(map(int,next(o)[1].split()))
    NIONS = sum(ntypes)
    if len(atoms) == 0: atoms = range(NIONS)
    ref_coords = np.zeros((len(atoms,3))
    old_coords = np.zeros((len(atoms,3))
    new_coords = np.zeros((len(atoms,3))
    displacement = np.zeros((len(atoms),3))

    # Reference frame 0
    next(o);
    for i in range(NIONS):
        if i in atoms:
            ref_coords[atoms.index(i)] = np.matmul(list(map(float,next(o)[1].split())),lat)
    old_coords = np.copy(ref_coords)

    # frame ~
    for idx, line in o:
        for i in range(step*(NIONS+1)): next(o);
        next(o);
        for i in range(NIONS):
            line = next(o)[1]
            if i in atoms:
                tmp = np.matmul(list(map(float,line.split())),lat)
                res = get_DIST(lat, old_coords[atoms.index(i)], tmp)
                displacement[atoms.index(i)] += res[1]
                new_coords[atoms.index(i)] = old_coords[atoms.index(i)] + res[1]

        old_coords = np.copy(new_coords)
        MSD.append(np.mean((displacement)**2))

    plt.plot(MSD)
    plt.tight_layout()
    plt.show()

def get_DIST(lat, pos_i, pos_j):
    #pos_X : Cartesian coordinates
    dist = 999
    for a in range(-1,2):
        for b in range(-1,2):
            for c in range(-1,2):
                tmp = pos_j + np.matmul(np.array([a,b,c]),lat)
                eval_dist = np.linalg.norm(pos_i-tmp)
                if eval_dist < dist:
                    dist = eval_dist
                    vec = tmp-pos_i
    return dist, vec

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
