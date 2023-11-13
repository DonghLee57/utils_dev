# python3.X
import argparse
import sys
import numpy as np
from ase.io import read, write
import ase.build

def flag_align(arg):
    if isinstance(arg,bool):
        return arg
    elif arg.lower() in ['c','C','T','t','1']:
        return True
    elif arg.lower() in ['b','B','F','f','0']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser(description='Build surface slab using ASE')
parser.add_argument('infile',type=str,help="Input POSCAR name")
parser.add_argument('outfile',type=str,help="Output POSCAR name")
index = parser.add_argument('index',type=str,action='store',help="Miller index as string ex) '1 1 1'")
parser.add_argument('n',type=int,help="The number of layers")
parser.add_argument('-v',required=False,default=7.5,type=float,help="Thickness of vacuum slab")
parser.add_argument('-a','-align',type=flag_align,default='c',required=False, help="Slab aligned to center/bottm")

args = parser.parse_args()

h,k,l = list(map(int,args.index.split()))

poscar = read(args.infile,format='vasp')
surf = ase.build.surface(poscar, (h,k,l), args.n, vacuum = args.v)
surf = surf[surf.numbers.argsort()]
if not args.a:
    minZ = np.min(surf.positions[:,2])
    t_vec = [0,0,-minZ]
    surf.translate(t_vec)
else:
    maxZ = np.max(surf.positions[:,2])
    t_vec = [0,0,(surf.cell[2][2]-maxZ)]
    surf.translate(t_vec)
write(args.outfile, images=surf, format='vasp')
