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
    flag = 0
    r_crt = 1.6
    rho = 2.2
    num_conc = 1E21 # atoms/cm^3

    if flag == 0:
        types = ['Si','O']
        atoms = [50, 100]
        symbols = []
        for i in range(len(types)):
            for j in range(atoms[i]):
                symbols += [types[i]]
        my_str = randomSTR(rho, symbols, atoms, r_crt)

    elif flag == 1:
        base = read(sys.argv[1],format='vasp')
        V = np.linalg.det(base.cell)
        num_add = 0
        while True:
            if num_add/V*1E24 < num_conc:
                num_add += 1
                print(num_add/V*1E24, num_conc)
            else:
                num_add -= 1
                break
        print(num_add)
        types = ['H']
        atoms = [12]
        symbols = base.get_chemical_symbols()
        for i in range(len(types)):
            for j in range(atoms[i]):
                symbols+=[types[i]]
        my_str = insertSTR(base, rho, symbols, atoms, r_crt)

#    write('POSCAR_generated.vasp',images=my_str,format='vasp')

#--------------------------------------------------------------------------------------
def randomSTR(rho,symbols,atoms,r_crt):
    IterMAX = 1000
    TOT_MASS = 0
    for i in symbols:
        TOT_MASS += atomic_masses[atomic_numbers[i]]
    lat = (TOT_MASS/rho/Na*10**24)**(1./3.)
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

def insertSTR(Obj,rho,symbols,atoms,r_crt):
    IterMAX = 1000
    TOT_MASS = 0
    for i in symbols:
        TOT_MASS += atomic_masses[atomic_numbers[i]]
    lat = (TOT_MASS/rho/Na*10**24)**(1./3.)
    lattice = np.zeros((3,3))
    for i in range(3): lattice[i][i] = lat
    Obj.cell = lattice
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

def rho2num(rho):
  
#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
