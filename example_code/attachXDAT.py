# python3.x
import sys
import ase.io

xdat1=ase.io.read(sys.argv[1],index=':',format='vasp-xdatcar')
xdat2=ase.io.read(sys.argv[2],index=':',format='vasp-xdatcar')

XDAT = xdat1 + xdat2
ase.io.write('XDATCAR_attached',images=XDAT,format='vasp-xdatcar')
