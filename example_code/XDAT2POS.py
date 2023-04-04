# python 3.x
import sys
import ase.io

XDAT = ase.io.read(sys.argv[1],index=':',format='vasp-xdatcar')
if len(sys.argv[2:]) > 1:
    [start,end,step] = list(map(int,sys.argv[2:]))
    if end == -1: 
        end = len(XDAT)
     for i in range(start,end+1,step):
         ase.io.write(f'POSCAR_{i}',images=XDAT[i-1],format='vasp')
       
elif len(sys.argv[2:]) == 1:
    fr = int(sys.argv[2])
    ase.io.write(f'POSCAR_{fr}',images=XDAT[fr-1],format='vasp')
