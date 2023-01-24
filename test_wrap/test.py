import numpy as np
import ase.io
from ase import Atoms
from ase import neighborlist

slab = ase.io.read(gpath + 'TiN_slab_mol_1.poscar',index=':',format='vasp')
lat = slab[0].cell
ilat = np.linalg.inv(lat)
types = slab[0].get_chemical_symbols()

crt_h = 0.45
candidates = []
for idx, pos in enumerate(slab[0].positions):
  if np.matmul(pos,ilat)[2] > crt_h:
    if types[idx] == 'N':
      candidates.append(idx)

print(candidates)
for i in enumerate(candidates):
  cut = [1.]*len(slab[0].positions)
  nei =  neighborlist.NeighborList(cut, self_interaction=False,bothways=True)
  nei.update(slab[0])
  nei_idx = nei.get_neighbors(i[1])[0]
  if len(nei_idx) > 5: candidates.remove(i[1])
  print(len(nei_idx))
print(candidates)


if len(candidates) > 0 :
  S = np.random.randint(0,len(candidates))
  cut = [1.]*len(slab[0].positions)
  nei =  neighborlist.NeighborList(cut, self_interaction=False,bothways=True)
  nei.update(slab[0])
  nei_idx = nei.get_neighbors(candidates[S])[0]
  vec = np.array([0,0,0])
  vec2 = np.array([0,0,0])
  for idx, i in enumerate(nei_idx):
    vec = vec + slab[0].get_distance(candidates[S], i, mic=True, vector=True)
  vec = vec/np.linalg.norm(vec)#*(-1)*2.2
  print(candidates[S], nei_idx)


mol = ase.io.read(gpath + 'TiCl4_mol.poscar',index=':',format='vasp')
mol[0].positions = mol[0].positions - np.mean(mol[0].positions,axis=0)
mol_axis = mol[0].positions[1] - mol[0].positions[0]
theta = vec_ang(vec, mol_axis)


#print(mol_axis)
#print(theta*180/np.pi)
### how to set universal molecule axis??

rot_axis = np.cross(vec, mol_axis)

### theta, coordinates sign...
mol[0].positions = rot_coords(mol[0].positions, theta, rot_axis) + slab[0].positions[candidates[S]] + vec*(-1)*2.2
mol[0].pop(1)
tmp_slab = slab.copy()
tmp_slab[0] = tmp_slab[0] + mol[0]


# check overlap
tag = True
#overlapping = slab[0].get_distances(mol[0], indices=range(len(slab[0].positions)), mic=True)
for idx, pos in enumerate(tmp_slab[0].positions):
  if idx < len(tmp_slab[0].positions) - len(mol[0].positions):
    overlapping = tmp_slab[0].get_distances(idx, indices=[-4,-3,-2,-1], mic=True)
    print(idx, min(overlapping))
    if min(overlapping) < 1.5:
      tag = False

if tag:
  ase.io.write(gpath+'TiN_slab_mol_tmp.poscar', images=tmp_slab, format='vasp') 
