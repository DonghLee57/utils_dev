# python 3.X
import sys
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.constraints import FixAtoms
from ase.data import atomic_numbers, atomic_names, atomic_masses
import scipy.constants as CONST
Na = CONST.N_A

def main():
    # 0: Generatating a random structure with the given condition (density, minimum distance)
    # 1: Inserting atoms based on an input structure
    flag = 1
    r_crt = 1.6
    rho = 2.2 # g/cm^3
    num_conc = 1E21  # atoms/cm^3

    if flag == 0:
        types = ['Si','O']
        atoms = [50, 100]
        symbols = []
        for i in range(len(types)):
            for j in range(atoms[i]):
                symbols += [types[i]]
        my_str = randomSTR(symbols, atoms, r_crt, rho=rho)
        #my_str = randomSTR(symbols, atoms, r_crt, numdensity=num_conc)

    elif flag == 1:
        base = read(sys.argv[1],format='vasp')
        num_add = get_num_addatoms(base, num_conc)
        print(num_add)
        types = ['H']
        atoms = [num_add[0]]
        symbols = base.get_chemical_symbols()
        for i in range(len(types)):
            for j in range(atoms[i]):
                symbols+=[types[i]]
        my_str = insertSTR(base, symbols, atoms, r_crt)

    write('POSCAR_generated.vasp',images=my_str,format='vasp')

#--------------------------------------------------------------------------------------
def randomSTR(symbols, atoms, r_crt, rho=None, numdensity=None):
    IterMAX = 1000
    if rho != None:
        TOT_MASS = 0
        for i in symbols:
            TOT_MASS += atomic_masses[atomic_numbers[i]]
        lat = (TOT_MASS/rho/Na*10**24)**(1./3.)
    if numdensity != None:
        lat = (len(symbols)/numdensity*1E24)**(1./3.)
    lattice = np.zeros((3,3))
    for i in range(3): lattice[i][i] = lat
    Obj = None
    for i in range(len(symbols)):
        trials = 0
        if Obj == None:
            Obj = Atoms(symbols[i],positions=[[0,0,0]],cell=lattice, pbc=[True,True,True])
            constraint = [len(Obj.positions)]
        else:
            while True and IterMAX > trials:
                trials += 1
                Obj.append(symbols[i])
                Obj.positions[-1] = np.matmul(np.random.random(3),Obj.cell)
                min_dist = min(Obj.get_distances(-1, indices=list(range(len(Obj.positions)-1)),mic=True))
                if r_crt > min_dist:
                   del Obj[-1]
                else: break
    Obj.set_constraint(FixAtoms(mask=constraint))
    return Obj

def insertSTR(Obj, symbols, atoms, r_crt):
    IterMAX = 1000
    for i in range(sum(atoms)):
        trials = 0
        while True and IterMAX > trials:
             trials += 1
             Obj.append(symbols[-sum(atoms)+i])
             Obj.positions[-1] = np.matmul(np.random.random(3),Obj.cell)
             min_dist = min(Obj.get_distances(-1, indices=list(range(len(Obj.positions)-1)),mic=True))
             if r_crt > min_dist:
                 del Obj[-1]
             else: break
    return Obj

def get_num_addatoms(Obj, NUMC):
    V = np.linalg.det(Obj.cell)
    ADD = 0
    ARR = []
    while True:
        diff =  ADD/V*1E24 - NUMC
        if diff < 0 :
            ADD += 1
            ARR.append(np.fabs(diff))
        else: break
    ARR = np.array(ARR)
    return [np.argmin(ARR),np.min(ARR)]
#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
