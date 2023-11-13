from ase.io import read, write
import ase.build
import numpy as np
# 'C' & 'H' represent cubic diamond and hexagonal diamond stacking, repectively.
# Recommend that the condition ((nC % 3 == 0) and (nH % 2 == 0)) should be satisfactory.
seq = ['C','C','H','H']

# Pre-defined Parameters
nC, nH = 1, 0
UNIT = 'POSCAR_diamond_unit'
layer_z = 3.1573715829248457 # angstrom
shift_vec = np.array([1/6,-1/6,0])

def flip(unit):
    t_vec = np.mean(unit.positions[:,2])
    unit.translate([0,0,-t_vec])
    unit.positions = unit.positions*np.array([1,1,-1])
    unit.translate([0,0,t_vec])
    return unit

def mk_unit_stack(unitcell):
    surf = ase.build.sruface(unitcell, (1,1,1), 1, vacuum = 7.5)
    surf = surf[surf.numbers.argsort()]
    minZ = np.min(surf.positions[:,2])
    surf.translate([0,0,-minZ])
    return surf
# 
unit = mk_unit_stack(read(UNIT))
lat = unit.cell
stack_vec = shift_vec@lat
stack_vec[2] = layer_z
unit.cell[2][2] = (len(seq)+1)*layer_z

NEW = unit.copy()
for idx, item in enumerate(seq):
    if item == 'C':
        nC += 1
        unit.translate(stack_vec)
        NEW += unit.copy()
    elif item == 'H':
        nH += 1
        unit = flip(unit)
        unit.translate([0,0,layer_z])
        stack_vec = stack_vec*np.array([-1,-1,1])
        NEW += unit.copy()
print(f"The number of cubic diamond layer is {nC:d}.")
print(f"The number of hexagonal diamond layer is {nH:d}.")
if (nC % 3 == 0) and (nH % 2 == 0):
    print(f"The generated stack is well matched.")
else:
    print(f"The generated stack may be mismatched.")
NEW.wrap()
write('POSCAR_stack.vasp',images=NEW,format='vasp')
