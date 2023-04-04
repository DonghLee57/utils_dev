import numpy as np
import sys
import ase.io
from ase import Atoms

base = ase.io.read(sys.argv[1],format='vasp')
NIONS = len(base.positions)
COM = base.get_center_of_mass()
MID = np.sum(base.cell,axis=0)/2
symbols = base.get_chemical_symbols()
base.positions = base.positions - COM + MID
base = base[base.numbers.argsort()]
ase.io.write(f'POSCAR_center', images=base, format='vasp')
