import sys
import numpy as np
from ase import Atoms
from ase.io import read

r_max= 6
dr = 0.02

base = read(sys.argv[1],format='vasp')
volume = np.fabs(np.linalg.det(base.cell))
positions = base.positions.copy()
symbols = np.array(base.get_chemical_symbols())
types = list(set(symbols))

combinations = []
for idx, i in enumerate(types):
    for jdx, j in enumerate(types):
        if idx <= jdx:
            combinations.append((i,j))

bins = np.arange(0.01, r_max, dr)
for idx, pair in enumerate(combinations):
    i = np.where(symbols == pair[0])[0]
    j = np.where(symbols == pair[1])[0]
    dist = []
    for a in i:
        for b in j:
            dist.append(base.get_distance(a,b,mic=True))
    res = np.histogram(dist,bins=bins)
    plt.plot(res[1][:-1],volume*res[0]/(len(i)*len(j)*4*np.pi*res[1][:-1]**2*dr),label=f'{pair[0]}-{pair[1]}')
plt.legend()
plt.show()
