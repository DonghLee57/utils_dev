import sys
import numpy as np
from ase.io import read, write
from ase import Atoms
from scipy.constants as CONST

#
IDX4REP = np.array([0,1,2,3])
CIDX_REF= -1
IDX4RM  = 0
CIDX_MOL= -1

def main():
    temp_anneal, time_anneal = 300, 10 # K, ps

    my = SIMULATOR()

    # initialization
    my.substrate('POSCAR_substrate.vasp')
    my.precursor('POSCAR_precursor.vasp')

    my.passivation('H', 3, prob=1.0)
    my.attach()
    my.optimize()

    my.determine()
    my.anneal(temp_anneal, t_anneal)

#
class SIMULATOR(T:float):
    ERRLOG = open('ERRLOG','w')
    q = CONST.e
    kb_eV = CONST.k/q
    def __init__(self) -> None:
        self.substrate = None
        self.precursor = None
        self.temperature = T

    def substrate(self, POSCAR:str) -> None:
        if self.substrate == None: self.substrate = read(POSCAR)
        else: ERRLOG.write('Your substrate is already assigned.')

    def precursor(self, POSCAR:str, center:int, neighbors) -> None:
        if self.precursor == None:
            self.precursor = read(POSCAR)
	    self.pre_center = center
	    self.pre_nei = neighbors
        else: ERRLOG.write('Your precursor is already assigned.')

    def passivate(self, elem:str, depth:float, prob:float=1.0):
        if self.substrate != None:
            pass
        else: ERRLOG.write('Your substrate is not assigned.')
            
    def optimize(self) -> None:
        if self.substrate != None:
            pass
        else: ERRLOG.write('Your substrate is not assigned.')

    def attach(self):
        if self.substrate != None:
            pass
        else: ERRLOG.write('Your substrate is not assigned.')
        if self.precursor != None:
            pass
        else: ERRLOG.write('Your precursor is not assigned.')

    def determine(self, temperature) -> None:
        oldE = calculate_energy(self.old_structure)
        newE = calculate_energy(self.new_structure)
        delta = newE - oldE
        if delta < 0:
            self.old_structure = self.new_structure
        else:
            prob = np.exp( -( delta ) / ( kb_eV * self.temperature ) )
            if np.random.rand() <  prob:
                self.old_structure = self.new_structure

    def anneal(self, temperature, time) -> None:
        pass
#
"""
def chemisorption():
    REF_STR = read(sys.argv[1],format='vasp')

    # Replace the assinged bonds
    NEW = Atoms('X',cell=REF_STR.cell, pbc=[1,1,1])
    NEW += REF_STR[CIDX_REF]
    for bidx, bond in enumerate(IDX4REP):
        SUB_MOL = read(sys.argv[2],format='vasp')
        SUB_MOL.translate(-SUB_MOL.get_center_of_mass())
        AX_m = SUB_MOL.positions[CIDX_MOL] - SUB_MOL.positions[IDX4RM]
        del SUB_MOL[IDX4RM]
        AX_b = REF_STR.positions[CIDX_REF] - REF_STR.positions[bond]
        AX_r = np.cross(AX_m, AX_b)
        theta = np.arccos(np.round(np.dot(AX_m,AX_b)/np.linalg.norm(AX_m)/np.linalg.norm(AX_b),4))
        if theta < 1E-4:
            SUB_MOL.positions = (-1)*(SUB_MOL.positions)
        else:
            SUB_MOL.positions = (-1)*rot_coords(SUB_MOL.positions, -theta, AX_r)
        SUB_MOL.translate(REF_STR.positions[bond])
        NEW += SUB_MOL.copy()

    # Write the output       
    del NEW[0]
    NEW = NEW[NEW.numbers.argsort()]
    write('POSCAR_replace.vasp',images=NEW, format='vasp')
"""
# Funtions
def rot_coords(coords, theta, axis):
    # theta: radian
    # axis: arbitarary axis to rotate
    ux, uy, uz = axis/np.linalg.norm(axis)
    R = np.array([[np.cos(theta)+ux**2*(1-np.cos(theta)),       ux*uy*(1-np.cos(theta))-uz*np.sin(theta),    ux*uz*(1-np.cos(theta))+uy*np.sin(theta)],\
                  [uy*ux*(1-np.cos(theta))+uz*np.sin(theta),    np.cos(theta)+uy**2*(1-np.cos(theta)),       uy*uz*(1-np.cos(theta))-ux*np.sin(theta)],\
                  [uz*ux*(1-np.cos(theta))-uy*np.sin(theta),    uz*uy*(1-np.cos(theta))+ux*np.sin(theta),    np.cos(theta)+uz**2*(1-np.cos(theta))]])
    return np.matmul(coords,R)

def vec_ang(v1, v2, exp='rad'):
    # exp: 'rad' or 'deg'
    angle = np.arccos(np.round(np.fabs(np.dot(v1,v2))/np.linalg.norm(v1)/np.linalg.norm(v2),4))
    if exp == 'rad':
        return angle
    elif exp == 'deg':
        return angle*180/np.pi

def fibonacci_sphere(samples:int=1, radius:float=1):
    """
    Generates points on a sphere using Fibonacci grid.

    :param samples: Number of points to generate.
    :param radius: Radius of the sphere.
    :return: Array of points on the sphere.
    """
    points = np.zeros((samples, 3))
    phi = np.pi * (3. - np.sqrt(5.))  # Golden angle

    for i in range(samples):
        y =  1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        R = np.sqrt( 1 - y**2)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * R
        z = np.sin(theta) * R

        points[i] = np.array([x, y, z])*radius
    return np.array(points)


def find_optimal_passivation_direction(IDX, Obj, radius=1.0, samples=1, min_distance=0.5):
    """
    Find the optimal direction for passivation on an atom.
    :param atom: Atom to passivate.
    :param atoms: ASE Atoms object representing the structure.
    :param radius: Radius for the Fibonacci grid.
    :param samples: Number of points to sample on the sphere.
    :param min_distance: Minimum allowed distance from other atoms.
    :return: Optimal direction vector for passivation, None if no valid direction found.
    """
    fibonacci_points = fibonacci_sphere(samples, radius)
    optimal_direction = np.zeros(3)

    for ii, point in enumerate(fibonacci_points):
        test_position = Obj.positions[IDX] + point
        test_atom = Atoms(['X'], positions=[test_position])
        test_structure = Obj + test_atom

        distances = test_structure.get_distances(-1, range(test_structure.get_global_number_of_atoms()-1), mic=True)
        if len(np.where( distances < min_distance )[0]) == 0:
            optimal_direction += test_structure.get_distance(IDX, -1, mic=True, vector=True)
            
    return optimal_direction/np.linalg.norm(optimal_direction)*2

#
if __name__ == "__main__":
    main()
