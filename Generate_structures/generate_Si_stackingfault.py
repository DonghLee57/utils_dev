from ase.io import read, write
from ase import Atoms
import numpy as np
# 'C' & 'H' represent cubic diamond and hexagonal diamond stacking, repectively.
# Recommend that the condition ((nC % 3 == 0) and (nH % 2 == 0)) should be satisfactory.
seq = ['C','C','H','H']
#seq = ['H','H','H']
#seq = ['C','C']

# Pre-defined Parameters
nC, nH = 1, 0
layer_z = 3.1573715829248457 # angstrom
shift_vec = np.array([1/6,-1/6,0])

def flip(unit):
    t_vec = np.mean(unit.positions[:,2])
    unit.translate([0,0,-t_vec])
    unit.positions = unit.positions*np.array([1,1,-1])
    unit.translate([0,0,t_vec])
    return unit

def mk_unit_stack(lat_a, lat_c):
    lat_A = lat_a*np.array([[1   ,           0],
                            [-0.5,np.sqrt(3)/2]])
    lat = np.block([[lat_A,np.zeros((2,1))],
                    [np.zeros((1,2)),lat_c]])
    for n in range(8):
        if n == 0: stack = Atoms('Si')
        else:      stack.append('Si')
    stack.cell = lat
    stack.positions[:4] = np.array([[0.0, 0.0, 0],
                                    [0.5, 0.0, 0],
                                    [0.0, 0.5, 0],
                                    [0.5, 0.5, 0]])
    stack.positions[4:] = stack.positions[:4].copy() + np.array([1/6, -1/6, 0])
    stack.positions = stack.positions@lat
    stack.positions[4:] += np.array([0,0,0.785731])
    stack.wrap()
    return stack
# 
unit = mk_unit_stack(7.6985577253335729, 10)
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
