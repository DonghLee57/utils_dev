import sys
import json
import numpy as np
from ase import Atoms
from ase.io import  write
from ase.data import chemical_symbols

with open(sys.argv[1],'r') as jo:
    tmp = json.load(jo)

lattice = float(sys.argv[2])

out = Atoms(pbc=[1,1,1])
out.cell = np.eye(3)*lattice
elem = tmp['PC_Compounds'][0]['atoms']['element']
positions = np.array([tmp['PC_Compounds'][0]['coords'][0]['conformers'][0]['x'],
                      tmp['PC_Compounds'][0]['coords'][0]['conformers'][0]['y'],
                      tmp['PC_Compounds'][0]['coords'][0]['conformers'][0]['z']]).T

for i in range(len(elem)):
    out.append(chemical_symbols[elem[i]])
    out.positions[-1] = positions[i]

write('POSCAR_converted', images=out, format='vasp')
