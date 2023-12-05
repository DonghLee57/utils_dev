import numpy as np
from ase.io import read, write
import ase.build.supercells

L = 100 # Angstrom
my_symbol = ['A','B']

CORE = read('POSCAR_unit',format='vasp')
nL = L//CORE.cell[0][0] + 1
nZ = 50//CORE.cell[2][2] + 1
P = np.array([[nL, 0, 0],
              [0, nL, 0],
              [0, 0, nZ]])
sCORE = ase.build.supercells.make_supercell(CORE, P)
sCORE = sCORE[sCORE.numbers.argsort()]
COM = sCORE.get_center_of_mass()
COM[2] = 0
sCORE.translate(-COM)

MID = np.sum(sCORE.cell[:2,:2],axis=0)/2
ID = np.where( np.linalg.norm((sCORE.positions[:,:2]),axis=1) > L//2)[0]
del sCORE[ID]
sCORE = sCORE[sCORE.numbers.argsort()]
#write('POSCAR_test.vasp',images=sCORE,format='vasp')

SHELL = read('res.lammps', style='atomic', format='lammps-data')
sSHELL = ase.build.supercells.make_supercell(SHELL, np.eye(3)*3)
COM = sSHELL.get_center_of_mass()
COM[2] = 0
sSHELL.translate(-COM)
sSHELL = sSHELL[sSHELL.numbers.argsort()]
MID = np.sum(sSHELL.cell[:2,:2],axis=0)/2
ID = np.where( np.linalg.norm((sSHELL.positions[:,:2]),axis=1) < L//2+1)[0]
del sSHELL[ID]
sSHELL = sSHELL[sSHELL.numbers.argsort()]
sym = sSHELL.get_chemical_symbols()
for idx, item in enumerate(sym):
    if   item == 'H' : sym[idx] = my_symbol[0]
    elif item == 'He': sym[idx] = my_symbol[1]
sSHELL.set_chemical_symbols(sym)
#write('POSCAR_test.vasp',images=sSHELL,format='vasp')

NEW = sSHELL + sCORE
write('POSCAR_core-shell.vasp',images=NEW,format='vasp')

