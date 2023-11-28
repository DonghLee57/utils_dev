import numpy as np
from ase.io import read, write
from ase import Atoms

TRJ = [Atoms('X',pbc=(1,1,1))]
TRJ[0].append('X')
TRJ[0].cell = np.eye(3)*10
TRJ[0].cell[2][0] = 1
uTRJ = [Atoms('X')]
uTRJ[0].append('X')
uTRJ[0].cell = np.eye(3)*10
uTRJ[0].cell[2][0] = 1

N = 50
history = np.zeros((N,2,3))
for n in range(1,N):
    uObj  = Atoms('X')
    uObj.append('X')
    uObj.cell = np.eye(3)*10
    uObj.cell[2][0] = 1
    Obj  = Atoms('X',pbc=(1,1,1))
    Obj.append('X')
    Obj.cell = np.eye(3)*10
    Obj.cell[2][0] = 1
    for a in range(2):
        xyz = np.random.random(3)
        history[n][a] = history[n-1][a] + xyz
        uObj.positions[a] = history[n][a]
        Obj.positions[a] = history[n][a]
    uTRJ.append(uObj)
    TRJ.append(Obj)

write('wrap.vasp',images=TRJ,format='vasp-xdatcar')
write('unwrap.vasp',images=uTRJ,format='vasp-xdatcar')
