import numpy as np
import sys, os
from ase.data import atomic_masses, atomic_numbers
import scipy.constants as CONST

energy = 10
mol = ['H','F']
coeff_mol= {'H': 1,
            'F': 1}

mass = 0
for m, sym in enumerate(mol):
  mass += atomic_masses[atomic_numbers[sym]]*coeff_mol[sym]

#Conversion
# g/mol >> kg
print(f"{mass:.2f} (g/mol)")
mass   *= 1E-3/(CONST.N_A)

# eV >> J
print(f"{energy:.2f} (eV)")
energy *= CONST.e

# m/s >> Angstrom/ps
velocity = np.sqrt(2*energy/mass)
print(f"{velocity:.2f} (m/s)")

velocity *= 1E-10/1E-12
print(f"{velocity:.2f} (Angs/ps)")
