import sys
import numpy as np
from ase.io import read
# For building a neighbor list
Rc  = 2.5

REF = read(sys.argv[1],format='vasp')
NEW = read(sys.argv[2],format='vasp')
NIONS = REF.get_global_number_of_atoms()

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
