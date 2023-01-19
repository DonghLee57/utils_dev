myLAMMPS = '/path/lammps/python'
import sys
sys.path.append(myLAMMPS)

import lammps
import ase.io

poscar = ase.io.read(sys.argv[1],index=':',format='vasp')
ase.io.write(sys.argv[2], images=poscar, format='lammps-data')

lmp = lammps.lammps()
lmp.file('script.in')
