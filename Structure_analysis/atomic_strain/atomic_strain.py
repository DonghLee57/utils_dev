import sys
import numpy as np
from ase.io import read
# For building a neighbor list
Rc  = 2.5

def match_order(REF, NEW):
    if REF.get_global_number_of_atoms() == NEW.get_global_number_of_atoms():
        NIONS = REF.get_global_number_of_atoms()
    else:
        print("!!! The size of two systems are different. !!!")
        return 1
    scaled_REF = REF.get_scaled_positions()
    scaled_NEW = NEW.get_scaled_positions()
    order = []
    for i in range(NIONS):
        for x in range(-1,2):
            for y in range(-1,2):
                for z in range(-1,2):
                    trans_vec = np.array([x,y,z])
                    dist = np.linalg.norm(scaled_NEW - (scaled_REF[i]+trans_vec),axis=1)
                    if np.fabs(np.min(dist)) < 0.01:
                        DIST = dist.copy()
        #print(f"{i:<4d} {np.fabs(np.min(DIST)):.4f} {np.min(DIST):.4f}")
        order.append(np.argmin(np.linalg.norm(scaled_NEW - scaled_REF[i],axis=1)))
    NEW = NEW[np.array(order)]
    return NEW

REF = read(sys.argv[1],format='vasp')
NEW = read(sys.argv[2],format='vasp')
NIONS = REF.get_global_number_of_atoms()
REF = match_order(NEW, REF)

All_E = np.empty((NIONS, 3, 3))
for idx in range(NIONS):
    X0 = REF.get_distances(idx, range(NIONS), mic=True,vector=True)
    dist = np.linalg.norm(X0,axis=1)
    IDs = np.where(dist < Rc)[0]
    IDs = IDs[np.where(IDs != idx)[0]]
    X0 = X0[IDs]

    X = NEW.get_distances(idx, range(NIONS), mic=True,vector=True)
    dist = np.linalg.norm(X,axis=1)
    IDs = np.where(dist < Rc)[0]
    IDs = IDs[np.where(IDs != idx)[0]]
    X = X[IDs]

    V = np.zeros((3,3))
    for jdx, item in enumerate(X0):
        V += np.matmul(np.array([X0[jdx]]).T,np.array([X0[jdx]]))
    W = np.zeros((3,3))
    for jdx, item in enumerate(X0):
        W += np.matmul(np.array([X0[jdx]]).T,np.array([X[jdx]]))
    J = np.matmul(np.linalg.inv(V),W)
    E = 0.5*(np.matmul(J,J.T)-np.eye(3))
    All_E[idx] = E

# OUTPUT: xyz format with strain
NEW.new_array('eps_xx', All_E[:,0,0],dtype=np.float64)
NEW.new_array('eps_yy', All_E[:,1,1],dtype=np.float64)
NEW.new_array('eps_zz', All_E[:,2,2],dtype=np.float64)
write('Analysis.xyz', NEW, format='extxyz')
