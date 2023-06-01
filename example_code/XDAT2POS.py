# python3.X
import sys
import ase.io

XDAT = ase.io.read(sys.argv[1],index=':',format='vasp-xdatcar')
if len(sys.argv[2:]) > 1:
    [start,end,step] = list(map(int,sys.argv[2:]))
    if end == -1: 
        end = len(XDAT)
    for i in range(start+step-1,end,step):
        ase.io.write(f'POSCAR_{i+1}',images=XDAT[i],format='vasp')

"""
elif len(sys.argv[2:]) == 1:
    # fr : 1 ~ max_image
    fr = int(sys.argv[2])
    if fr == -1: 
        fr = len(XDAT)
    ase.io.write(f'POSCAR_{fr}',images=XDAT[fr-1],format='vasp')
"""
