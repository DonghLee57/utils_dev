import sys
import numpy as np
from ase.io import read, write
from ase import Atoms
import scipy.constants as CONST

def main():
    temp_anneal, time_anneal = 300, 10 # K, ps

    my = SIMULATOR(temp_anneal)

    # initialization
    my.set_substrate('POSCAR_substrate.vasp')
    my.set_precursor('POSCAR_precursor.vasp', -1, list(range(4)))

    psv = my.passivate('H', 1.5, 1.5)
    psv = psv[psv.numbers.argsort()]
    write("psv.vasp", images=psv, format='vasp')

    sym = np.array(psv.get_chemical_symbols())
    condition_a = sym != 'H'
    condition_b = psv.positions.T[2] > (np.max(psv.positions.T[2]) - 1.5)
    condition = condition_a*1 + condition_b
    ids = np.where( condition == 2 )[0]
    ids = np.random.choice(ids, size=len(ids)//4)
    for ii, item in enumerate(ids):
        my.attach(item)
    
    #my.optimize()

    #my.determine()
    #my.anneal(temp_anneal, t_anneal)

#
class SIMULATOR:
    ERRLOG = open('ERRLOG','w')
    q = CONST.e
    kb_eV = CONST.k/q
    def __init__(self, T:float) -> None:
        self.substrate = None
        self.precursor = None
        self.temperature = T
        self.r_overlap = 1.0

    def set_substrate(self, POSCAR:str) -> None:
        if self.substrate == None:
            self.substrate = read(POSCAR)
            self.old_structure = self.substrate.copy()
        else: ERRLOG.write('Your substrate is already assigned.')

    def set_precursor(self, POSCAR:str, center:int, neighbors) -> None:
        if self.precursor == None:
            self.precursor = read(POSCAR)
            self.precursor.translate(-self.precursor.get_center_of_mass())
            self.pre_center = center
            self.pre_nei = neighbors
        else: ERRLOG.write('Your precursor is already assigned.')

    def passivate(self, elem:str, length:float ,depth:float, prob:float=1.0):
        if self.old_structure != None:
            obj = self.old_structure
            ids = np.where( obj.positions.T[2] > np.max(obj.positions.T[2]) - depth )[0]
            for ii, item in enumerate(ids):
                dvec = find_opt_direction(obj, item, elem, length, 50, self.r_overlap)
                test_position = obj.positions[item] + dvec
                test_atom = Atoms([elem], positions=[test_position])
                test_structure = obj.copy() + test_atom
                distances = test_structure.get_distances(-1, range(test_structure.get_global_number_of_atoms()-1), mic=True)
                if (len(np.where( distances < self.r_overlap )[0]) == 0) and np.random.rand() < prob:
                    obj = obj + test_atom
            self.old_structure = obj.copy()
            return obj
        else:
            ERRLOG.write('Your substrate is not assigned.')
            sys.exit()
   
    def optimize(self):
        if self.substrate != None:
            pass
        else: ERRLOG.write('Your substrate is not assigned.')

    # debugging
    def attach(self, idx:int, elem:str):
        if self.old_structure != None:
            if self.precursor != None:
                obj = self.old_structure.copy()
                distances = obj.get_distances(idx, range(obj.get_global_number_of_atoms()), mic=True)
                # check
                jdx = np.random.choice(np.where( distances < 1.5 )[0], size=1)
                ax_b = obj.positions[idx] - obj.positions[jdx]
                
                mol = self.precursor.copy()
                # modify...
                ax_m = mol.positions[self.pre_center] - mol.positions[self.pre_nei[0]]
                del mol[self.pre_nei[0]]
                
                ax_r = np.cross(ax_m, ax_b)
                theta = np.arccos( np.round( np.dot(ax_m,ax_b)/np.linalg.norm(ax_m)/np.linalg.norm(ax_b) , 4) )
                if theta < 1E-4:
                    mol.positions = (-1)*(mol.positions)
                else:
                    mol.positions = (-1)*rot_coords(mol.positions, -theta, AX_r)
                mol.translate(obj.positions[jdx])

                self.old_structure = obj.copy()
            else:
                ERRLOG.write('Your precursor is not assigned.')
                sys.exit()
        else:
            ERRLOG.write('Your substrate is not assigned.')
            sys.exit()

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

# Funtions
def rot_coords(coords, theta, axis):
    #Rotate coordinates with angle theta along the axis.
    # coords: coordinates of atoms
    # theta: radian
    # axis: axis to rotate
    ux, uy, uz = axis/np.linalg.norm(axis)
    R = np.array([[np.cos(theta)+ux**2*(1-np.cos(theta)),       ux*uy*(1-np.cos(theta))-uz*np.sin(theta),    ux*uz*(1-np.cos(theta))+uy*np.sin(theta)],\
                  [uy*ux*(1-np.cos(theta))+uz*np.sin(theta),    np.cos(theta)+uy**2*(1-np.cos(theta)),       uy*uz*(1-np.cos(theta))-ux*np.sin(theta)],\
                  [uz*ux*(1-np.cos(theta))-uy*np.sin(theta),    uz*uy*(1-np.cos(theta))+ux*np.sin(theta),    np.cos(theta)+uz**2*(1-np.cos(theta))]])
    return np.matmul(coords,R)

def vec_ang(v1, v2, exp='rad'):
    #Calculate the angle between two vectors, v1 and v2.
    # v1, v2:  same size of vector arrays.
    # exp: 'rad' or 'deg'
    angle = np.arccos(np.round(np.fabs(np.dot(v1,v2))/np.linalg.norm(v1)/np.linalg.norm(v2),4))
    if exp == 'rad':
        return angle
    elif exp == 'deg':
        return angle*180/np.pi

def find_opt_direction(obj, idx:int, elem:str, radius:float=2.0, npts:int=20, r_overlap:float=1.0):
    #Find the optimal direction for passivation on an atom.
    # param obj: ASE Atoms object representing the structure.
    # param idx: Atom to passivate.
    # param elem: element for passivation.
    # param radius: Radius for the Fibonacci grid.
    # param npts: Number of points on the sphere.
    # param r_overlap: Minimum allowed distance from other atoms.
    # return: Optimal direction vector for passivation, [0,0,0] if no valid direction found.
    fibpts = np.zeros((npts, 3))
    phi = np.pi * (3. - np.sqrt(5.))  # Golden angle
    for i in range(npts):
        y =  1 - (i / float(npts - 1)) * 2
        R = np.sqrt( 1 - y**2)
        theta = phi * i
        x = np.cos(theta) * R
        z = np.sin(theta) * R
        fibpts[i] = np.array([x, y, z]) * radius
    fibpts = np.array(fibpts)
    direction = np.zeros(3)
    for ii, point in enumerate(fibpts):
        test_position = obj.positions[idx] + point
        test_atom = Atoms([elem], positions=[test_position])
        test_structure = obj + test_atom
        distances = test_structure.get_distances(-1, range(test_structure.get_global_number_of_atoms()-1), mic=True)
        if len(np.where( distances < r_overlap )[0]) == 0:
            direction += test_structure.get_distance(idx, -1, mic=True, vector=True)
    return direction/np.linalg.norm(direction)*radius
 
#
if __name__ == "__main__":
    main()
