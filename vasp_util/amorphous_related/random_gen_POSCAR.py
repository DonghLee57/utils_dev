# python 3.X
import sys
import numpy as np
import ase
import ase.io
from ase import Atoms
from ase.constraints import FixAtoms
from ase.data import atomic_numbers, atomic_names, atomic_masses

def main():
    # 0: generate from scratch
    # 1: insert based on an input structure
    flag = 1
    r_crt = 1.5
    
    base = ase.io.read(sys.argv[1],format='vasp')
    print(f'Density(input): {get_density(base):.6f} g/cm3')
    rho = get_density(base)
    if flag == 0:
        # Generate a random structure from scratch
        types = ['Si']
        atoms = [100]
        symbols = []
        for i in range(len(types)):
            for j in range(atoms[i]):
                symbols+=[types[i]]
        amor = randomPOSCAR(rho, symbols, atoms, r_crt)
    
    elif flag == 1:
        # Insert random-position atoms based on a given structure
        types = ['H']
        atoms = [10]
        symbols = base.get_chemical_symbols()
        for i in range(len(types)):
            for j in range(atoms[i]):
                symbols+=[types[i]]
        amor = insertPOSCAR(base, rho, symbols, atoms, r_crt)
        
    ase.io.write('POSCAR_generated.vasp',images=amor,format='vasp')

#--------------------------------------------------------------------------------------

def get_density(atom_obj):
    Na = 6.022E+23
    symbols = atom_obj.get_chemical_symbols()
    TOT_MASS = 0
    for i in symbols:
        TOT_MASS += atomic_masses[atomic_numbers[i]]
    # Unit : g/cm^3
    return TOT_MASS/atom_obj.get_volume()/Na*10**24

def randomPOSCAR(rho,symbols,atoms,r_crt):
    IterMAX = 1000
    Na = 6.022E+23
    TOT_MASS = 0
    for i in symbols:
        TOT_MASS += atomic_masses[atomic_numbers[i]]
    lat = (TOT_MASS/rho/Na*10**24)**(1./3.)
    lattice = np.zeros((3,3))
    for i in range(3): lattice[i][i] = lat
    genPOSCAR = None
    for i in range(len(symbols)):
        trials = 0
        if genPOSCAR == None:
            genPOSCAR = Atoms(symbols[i],positions=[[0,0,0]],cell=lattice, pbc=[True,True,True])
            constraint = [len(genPOSCAR.positions)]
        else:
            while True and IterMAX > trials:
                trials += 1
                genPOSCAR.append(symbols[i])
                genPOSCAR.positions[-1] = np.matmul(np.random.random(3),genPOSCAR.cell)
                min_dist = min(genPOSCAR.get_distances(-1, indices=list(range(len(genPOSCAR.positions)-1)),mic=True))
                #print(len(genPOSCAR.positions),trials,min_dist)
                if r_crt > min_dist:
                   del genPOSCAR[-1]
                else: break
    genPOSCAR.set_constraint(FixAtoms(mask=constraint))
    return genPOSCAR

def insertPOSCAR(base,rho,symbols,atoms,r_crt):
    IterMAX = 1000
    Na = 6.022E+23
    TOT_MASS = 0
    for i in symbols:
        TOT_MASS += atomic_masses[atomic_numbers[i]]
    lat = (TOT_MASS/rho/Na*10**24)**(1./3.)
    lattice = np.zeros((3,3))
    for i in range(3): lattice[i][i] = lat
    base.cell = lattice
    for i in range(sum(atoms)):
        trials = 0
        while True and IterMAX > trials:
             trials += 1
             base.append(symbols[-sum(atoms)+i])
             base.positions[-1] = np.matmul(np.random.random(3),base.cell)
             min_dist = min(base.get_distances(-1, indices=list(range(len(base.positions)-1)),mic=True))
             if r_crt > min_dist:
                 del base[-1]
             else: break
    return base

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
