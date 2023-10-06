import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.eos import EquationOfState as eos
from scipy import constants
q = constants.e
k_eV = constants.Boltzmann/q

dirs = './'
V, DFT = [], []
for i in range(1,13):
    base = read(f'{dirs}/POSCAR_{i}',format='vasp')
    V.append(base.get_volume())
    DFT.append(read(f'{dirs}/OUTCAR_{i}',format='vasp-out').get_potential_energy())
NIONS = base.get_global_number_of_atoms()
dft = eos(V,DFT,eos='birchmurnaghan')
v0,e0,B = dft.fit()
print('DFT')
print(f'Vmin : {v0/NIONS:6.2f}')
print(f'Emin : {e0/NIONS:6.2f}')
print(f'Bulk modulus: {B*q/1E-21:6.2f} GPa')

fig, ax = plt.subplots(figsize=(5,5))
ax.tick_params(axis='both', which='major', labelsize=15)
ax.plot(dft.v/NIONS,dft.e/NIONS,label='DFT')
ax.scatter(np.array(V)/NIONS,np.array(DFT)/NIONS,marker='s')
ax.set_xlabel(r'Volume ($\mathrm{\AA}$/atom)',fontsize=15)
ax.set_ylabel('Energy (eV/atom)',fontsize=15)
ax.legend()

plt.tight_layout()
plt.show()
