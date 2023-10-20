import sys
import numpy as np
from ase.io import read, write
import ase.build.supercells
from ase import Atoms
from ase.data import atomic_numbers, atomic_names, atomic_masses

def main():
    unit = read(sys.argv[1],format='vasp')
    R = float(sys.argv[2])
    Nx = R//unit.cell[0][0] + 2
    P = np.array([[Nx, 0, 0],
                  [0, Nx, 0],
                  [0, 0, Nx]])
    supercell = ase.build.supercells.make_supercell(unit, P)
    MID = np.sum(supercell.cell,axis=1)/2
    IDs = np.where(np.linalg.norm(supercell.positions-MID,axis=1) > R)[0]
    del supercell[IDs]
    supercell.cell = np.eye(3)*(2*R+10)
    supercell = supercell[supercell.numbers.argsort()]
    MID = np.sum(supercell.cell,axis=1)/2
    supercell.translate(-supercell.get_center_of_mass()+MID)
    print(f"Maximum distance from the center: {np.max(np.linalg.norm(supercell.positions-MID,axis=1)):.2f} Angstrom")
    write(f'POSCAR_cluster_{R}.vasp',images=supercell,format='vasp')
    passive = passivation(supercell)
    if passive != 0:
        write(f'POSCAR_cluster_{R}_passive.vasp',images=passive,format='vasp')
    else:
        print("Input other size of R")
    return 1

#--------------------------------------------------------------------------------------
def passivation(Obj,ele='H'):
    RH = 1.5
    Rc = 2.6
    passive = Atoms('X')
    all_idx = list(range(len(Obj.positions)))
    for idx, item in enumerate(Obj.positions):
        vec = Obj.get_distances(idx, all_idx, mic=True,vector=True)
        dist = np.linalg.norm(vec,axis=1)
        IDs = np.where(dist <Rc)[0]
        IDs = IDs[np.where(IDs != idx)[0]]
        CN = len(IDs)
        if CN < 4:
            if CN == 2:
                for pdx, pos in enumerate(vec[IDs]):
                    passive.append(ele)
                    passive.positions[-1] = item - rot_coords(pos*RH/Rc, np.pi/2, -np.sum(vec[IDs],axis=0))
            elif CN == 3:
                p_vec = np.zeros(3)
                for pdx, pos in enumerate(vec[IDs]):
                    p_vec += pos*RH/Rc
                passive.append(ele)
                passive.positions[-1] = item - p_vec
            elif CN == 1:
                print('Si has one Si-Si bond!!!')
                return 0
    del passive[0]
    TOT = Obj + passive
    TOT = TOT[TOT.numbers.argsort()]
    return TOT

def rot_coords(coords, theta, axis):
    # theta: radian
    # axis: arbitarary axis to rotate
    ux, uy, uz = axis/np.linalg.norm(axis)
    R = np.array([[np.cos(theta)+ux**2*(1-np.cos(theta)),       ux*uy*(1-np.cos(theta))-uz*np.sin(theta),    ux*uz*(1-np.cos(theta))+uy*np.sin(theta)],\
                  [uy*ux*(1-np.cos(theta))+uz*np.sin(theta),    np.cos(theta)+uy**2*(1-np.cos(theta)),       uy*uz*(1-np.cos(theta))-ux*np.sin(theta)],\
                  [uz*ux*(1-np.cos(theta))-uy*np.sin(theta),    uz*uy*(1-np.cos(theta))+ux*np.sin(theta),    np.cos(theta)+uz**2*(1-np.cos(theta))]])
    return np.matmul(coords,R)
#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
