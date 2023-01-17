# @shell
# $ export ASE_LAMMPSRUN_COMMAND=/path/to/lmp_binary
# 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6719897/
import ase
from ase.calculators.lammpsrun import LAMMPS


lmp = lammps()
