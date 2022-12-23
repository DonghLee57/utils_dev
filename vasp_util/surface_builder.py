# python3.X
import argparse
import sys
import ase.io
import ase.build

parser = argparse.ArgumentParser(description='Build surface slab using ASE')
parser.add_argument('infile',type=str,help="Input POSCAR name")
parser.add_argument('outfile',type=str,help="Output POSCAR name")
index = parser.add_argument('index',type=str,action='store',help="Miller index as string ex) '1 1 1'")
parser.add_argument('n',type=int,help="The number of layers")
parser.add_argument('-v',required=False,default=7.5,type=float,help="Thickness of vacuum slab")

args = parser.parse_args()

print(args)

h,k,l = list(map(int,args.index.split()))

poscar = ase.io.read(args.infile,index=':',format='vasp')
surf111 = ase.build.surface(poscar[0], (h,k,l), args.n, vacuum = args.v)
ase.io.write(args.outfile, images=surf111, format='vasp')                                                        
