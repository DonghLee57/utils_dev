import sys
from ase.io import read, write

opt ={'atom_sort': True}

FILE = sys.argv[1]
[start, end, step] = list(map(int,sys.argv[2:5]))
#[start, end, step] = [1000, 5000, 10]

my_slice = slice(start-1+step, end-1+step, step)
XDAT = read(FILE ,index=my_slice, format='vasp-xdatcar')

for idx, item in enumerate(range(start+step, end+1, step)):
    if opt['atom_sort']:
        XDAT[idx] = XDAT[idx][XDAT[idx].numbers.argsort()]
    POSCAR = write(f'./structures/POSCAR_{item}',images=XDAT[idx],format='vasp')
