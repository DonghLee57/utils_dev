import numpy as np
from ase.io import read, write
from ase import Atoms

# functions
def relax(obj):
    # relax obj
    # obj: ASE Atoms object
    # return the relaxed structure
    pass

def mk_fibonacci_pts(center, npts, distance):
    # generate fibonacci lattices with the given arguments
    # npts: the number of points
    # distance: the distance from the center position
    # return list of coordinates
    pts = np.empty((npts,3))
    return pts

def check_overlap(obj, cutoff=1.0) -> bool:
    # calculate distances from the input fibonacci lattice
    # if lattice within cutoff -> False, else -> True 
    # obj: ASE Atoms object
    if 1:
        return True
    else:
        return False
        
def neb(initial, final):
    return Ea

def get_distance(a, b) -> float:
    dist = np.linalg.norm(b-a)
    return dist

# initial structure
tmp = Atoms('Si', positions=[[5,5,5]], pbc=[1,1,1], cell=np.eye(3)*10)
write('POSCAR',images=tmp, format='vasp')

# load structure data
structure = read('POSCAR', format='vasp')

# initial position
structure.append('H')
structure.positions[-1] = [6,6,6]

# Relax initial structure
relax(structure)

# make Fibonacci lattice
pts = mk_fibonacci_pts()

# check the Fibonacci lattice
res = np.zeros((len(pts),2))
ii=0
event = []
for idx, point in enumerate(pts):
    tmp = structure.copy()
    tmp.positions[-1] = point
    if check_overlap(tmp):
        relax_tmp = relax(tmp)
        event.append(relax_tmp)
        dist = get_distance(structure.positions[-1], tmp.positions[-1])
        Ea = neb(structure, tmp)
        res[ii] = [idx, Ea, dist]
        ii += 1

# kmc
kb_eV = 8E-15 # boltzmann
T = 300 # K
rate = np.zeros(ii+1)
for i in range(ii):    rate[i] = np.exp(-res[i][1]/(kb_eV*T))
Q = np.sum(rate)
up = np.random.random()
delta_time = np.log(1/up)/Q
t += delta_time
structure = event[index]
#go to check the Fibonacci lattice

"""
running_sum = 0
my_array = np.array([1,2,3,4,5])
threshold = 4
for index, entry in enumerate(my_array):
    running_sum += entry
    if running_sum > threshold:
        break
if running_sum < threshold:
    index = -1
print(index, np.cumsum(my_array))
# event[index]

print(np.log(1/0.001)/1E+13)
"""
