# python3.X

import sys
import numpy as np
import ase.io
import ase.build.supercells
from ase.constraints import FixAtoms

unit = ase.io.read(sys.argv[1],format='vasp')
P = np.array([[5, 0, 0],
              [0, 5, 0],
              [0, 0, 5]])
supercell = ase.build.supercells.make_supercell(unit, P)
supercell = supercell[supercell.numbers.argsort()]
supercell.set_constraint(FixAtoms([-1]))
ase.io.write('POSCAR-supercell', images=supercell, format='vasp')
