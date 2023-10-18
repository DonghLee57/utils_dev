# python3.X
import sys
import numpy as np
from ase.io import read, write
import ase.build.supercells

unit = read(sys.argv[1],format='vasp')
[xx,xy,xz] = list(map(int,sys.argv[2:]))
MAT = np.array([[xx, 0, 0],
                [0, xy, 0],
                [0, 0, xz]])
scell = ase.build.supercells.make_supercell(unit, MAT)
scell = scell[scell.numbers.argsort()]
write('POSCAR_supercell.vasp', images=scell, format='vasp')
