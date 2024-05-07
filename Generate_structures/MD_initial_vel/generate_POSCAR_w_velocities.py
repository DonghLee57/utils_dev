# python3.x
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import sys
import numpy as np
import scipy.constants as CONST
kb_eV = CONST.k/CONST.e

FILE = 'POSCAR_w_velocities.vasp'
poscar = read(sys.argv[1], format='vasp')
temperature = 300.0 #K

MaxwellBoltzmannDistribution(poscar, temperature_K=temperature*kb_eV)
write(FILE, images=poscar, format='vasp')

with open(FILE, 'a+') as o:
    o.write("C # units: angstrom/fs\n")
    for _, item in enumerate(poscar.get_velocities()):
        vx, vy, vz = item
        o.write(f"{vx:>14.6E} {vy:>14.6E} {vz:>14.6E}\n")
