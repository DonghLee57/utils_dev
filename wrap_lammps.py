# Before using LAMMPS in python,
# Go to lammps directory. Then,
# cd ./src/
# make (gcc/g++/mpi) mode=shlib
myLAMMPS = '/path/lammps/python'
import sys
sys.path.append(myLAMMPS)

import lammps
lmp = lammps.lammps()
