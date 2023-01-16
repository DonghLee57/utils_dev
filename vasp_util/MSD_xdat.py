import sys
import numpy as np

def main():
    atoms = []
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
    
    # Reference frame 0
    next(o);
    for i in range(NIONS):
        if i in atoms:
            ref_coords[atoms.index(i)] = list(map(float,next(o)[1].split()))
    old_coords = np.copy(ref_coords)

    # frame ~
    for idx, line in o:
        for i in range(step*(NIONS+1)): next(o);
        next(o);
        for i in range(NIONS):
            line = next(o)[1]
            if i in atoms:
                #res = get_DIST(lat, old_coords[n],new_coords[n])
                new_coords[atoms.index(i)] = list(map(float,line.split()))
                
        
#new_coords[n] = res[1]
#MSD = np.mean((new_coords - ref_coords)**2)


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
                    close_pos = np.copy(tmp)
    return dist, close_pos

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
