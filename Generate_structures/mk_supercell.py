# python3.X
import sys
import numpy as np
from ase.io import read, write
import ase.build.supercells
unit = read(sys.argv[1],format='vasp')
P = np.array([[5, 0, 0],
              [0, 5, 0],
              [0, 0, 5]])
supercell = ase.build.supercells.make_supercell(unit, P)
supercell = supercell[supercell.numbers.argsort()]
write('POSCAR-supercell', images=supercell, format='vasp')
