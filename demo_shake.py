import sys
import numpy as np
import ase.io
from ase import Atoms
from ase import neighborlist

def shake(atom_obj, delta, dev, crit_overlap):
  NIONS = len(atom_obj.positions)
  cut = [crit_overlap*0.5]*NIONS
  nei =  neighborlist.NeighborList(cut, skin =0, self_interaction=False,bothways=True)
  nei.update(atom_obj)
  for n in range(NIONS):
    perturb = np.random.randn(3)
    perturb = perturb/np.linalg.norm(perturb)*np.random.normal(delta,dev,1)
    
    nei_idx = nei.get_neighbors(n)[0]
    if len(nei_idx) < 1:
      atom_obj.positions[n] += perturb
    else:
      for i in range(len(nei_idx)):
        atom_obj.positions[n] -= atom_obj.get_distance(n, nei_idx[i], mic=True,vector=True)*0.1
  return atom_obj

unit = ase.io.read(sys.argv[1],index=':',format='vasp')
P = np.identity(3)*3
supercell = make_supercell(unit,P)
supercell = supercell[supercell.numbers.argsort()]

ref = supercell
NIONS = len(supercell.positions)
c = 0
for x in range(100):
  new = shake(supercell, 0.1, 0.01, 1)
  ase.io.write(save_file, images=new, append=True, format='extxyz')
  ref = new
 
