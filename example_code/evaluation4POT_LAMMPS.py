import numpy as np
import os,subprocess
import ase.io
import ase.build.supercells
from scipy import constants
q = constants.e

LAMMPS='~/lmp_mpi'
style ='atomic'

pot_type='tersoff'
pot = 'Si.tersoff'

def main():
    opt_str = relax('POSCAR-unit')
    EOS(opt_str)
    defects_calc(opt_str)
    surf_calc(opt_str)
    return 0
#--------------------------------------------------------------------------------------
def relax(struct):
    opt_str = 'Si_relax.txt'
    unit = ase.io.read(struct,format='vasp')
    ase.io.write('Si.txt',images=unit, format='lammps-data')
    with open(f'./script.txt','w') as o:
        o.write('units       metal\n')
        o.write(f'atom_style  {style}\n')
        o.write('boundary    p p p\n')
        o.write('box         tilt large\n')
        o.write('read_data   "Si.txt"\n')
        set_pot(o)
        o.write('mass 1 28.0855\n')
        o.write('fix         1 all box/relax iso 0.0 vmax 0.001 fixedpoint 0 0 0\n')
        o.write('minimize 1e-6 1e-10 10000 100000\n')
        o.write(f'write_data {opt_str}')
    target = subprocess.check_output([LAMMPS,'-in','script.txt'],universal_newlines=True)
    return opt_str

def EOS(struct,show=False):
    struct = ase.io.read(struct,style=style,format='lammps-data')
    x = np.arange(0.95,1.06,0.01)
    V, E = np.zeros(len(x)), np.zeros(len(x))
    for i in range(len(x)):
        name = 'str.dat'
        new = struct.copy()
        P = np.eye(3)*x[i]
        new.cell = np.matmul(struct.cell,P)
        new.positions = np.matmul(struct.positions,P)
        ase.io.write(name, images=new, format='lammps-data')
        V[i] = np.linalg.det(new.cell)
        E[i] = oneshot(new)
    from ase.eos import EquationOfState as eos_fit
    eos = eos_fit(V,E,eos='birchmurnaghan')
    V0, E0, B = eos.fit()
    print(f"Equilibrium volume: {V0:.2f} Angs^3")
    print(f"Minimum energy: {E0:.2f} eV")
    print(f"Bulk modulus: {B*q*1.0e21:.2f} GPa")
    if show:
        eos.plot(show=True)
    return 0

def defects_calc(struct, pot=pot):
    struct = ase.io.read(struct,style=style,format='lammps-data')
    P = np.eye(3)*2
    supercell = ase.build.supercells.make_supercell(struct, P)
    supercell = supercell[supercell.numbers.argsort()]

    E0 = oneshot(supercell)
    u = E0/supercell.get_global_number_of_atoms()
    V = supercell.copy(); del V[-1];
    E_V = minimize(V) - (E0 - u)
    os.system('mv res.dat Si_V.dat')

    Int_T = supercell.copy()
    id_list = np.array([38, 39, 42, 43])-1
    Int_T.append('H')
    Int_T.positions[-1] = np.mean(supercell.positions[id_list],axis=0)
    #ase.io.write('INT_T_Si',images=Int_T,format='lammps-data')
    E_T = minimize(Int_T) - (E0 + u)
    os.system('mv res.dat Si_T.dat')

    Int_H = supercell.copy()
    id_list = np.array([11, 39, 42, 43, 48, 57])-1
    Int_H.append('H')
    Int_H.positions[-1] = np.mean(supercell.positions[id_list],axis=0)
    #ase.io.write('INT_H_Si',images=Int_H,format='lammps-data')
    E_H = minimize(Int_H) - (E0 + u)
    os.system('mv res.dat Si_H.dat')

    print(f"------------Defect energy------------")
    print(f"Defect energy for V_Si : {E_V:.2f} eV")
    print(f"Defect energy for Si_T : {E_T:.2f} eV")
    print(f"Defect energy for Si_H : {E_H:.2f} eV")
    return 0        

def surf_calc(struct, pot=pot):
    struct = ase.io.read(struct,style=style,format='lammps-data')
    E0 = oneshot(struct)
    u = E0/struct.get_global_number_of_atoms()
    surf100 = ase.build.surface(struct, (1,0,0), 4, vacuum = 7.5)
    surf100 = surf100[surf100.numbers.argsort()]
    S100 = (minimize(surf100) - u*surf100.get_global_number_of_atoms())/(2*np.linalg.det(surf100.cell[:2,:2]))
    os.system('mv res.dat Si_100.dat')

    surf110 = ase.build.surface(struct, (1,1,0), 4, vacuum = 7.5)
    surf110 = surf110[surf110.numbers.argsort()]
    S110 = (minimize(surf110) - u*surf110.get_global_number_of_atoms())/(2*np.linalg.det(surf110.cell[:2,:2]))
    os.system('mv res.dat Si_110.dat')

    surf111 = ase.build.surface(struct, (1,1,1), 4, vacuum = 7.5)
    surf111 = surf111[surf111.numbers.argsort()]
    S111 = (minimize(surf111) - u*surf111.get_global_number_of_atoms())/(2*np.linalg.det(surf111.cell[:2,:2]))
    os.system('mv res.dat Si_111.dat')

    print(f"-----------Surface energy-----------")
    print(f"(100) Surface energy : {S100:.2f} eV")
    print(f"(110) Surface energy : {S110:.2f} eV")
    print(f"(111) Surface energy : {S111:.2f} eV")
    return 0

def oneshot(struct, pot=pot, pot_type=pot_type):
    name = 'str.dat'
    ase.io.write(name,images=struct,format='lammps-data')
    with open(f'./script.txt','w') as o:
        o.write('units       metal\n')
        o.write(f'atom_style  {style}\n')
        o.write('boundary    p p p\n')
        o.write('box         tilt large\n')
        o.write(f'read_data   {name}\n')
        set_pot(o)
        o.write('mass 1 28.0855\n')
        o.write('run 1\n')
    target = subprocess.check_output([LAMMPS,'-in','script.txt'],universal_newlines=True).split('\n')
    o=enumerate(target)
    for idx, line in o:
        tmp = line.split()
        if len(tmp) > 0:
            if tmp[0] =='Step':
                next(o)
                step, T, PE = list(map(float,next(o)[-1].split()[:3]))
    return PE    

def minimize(struct, pot=pot, pot_type=pot_type):
    name = 'str.dat'
    ase.io.write(name,images=struct,format='lammps-data')
    with open(f'./script.txt','w') as o:
        o.write('units       metal\n')
        o.write(f'atom_style  {style}\n')
        o.write('boundary    p p p\n')
        o.write('box         tilt large\n')
        o.write(f'read_data   {name}\n')
        set_pot(o)
        o.write('mass 1 28.0855\n')
        o.write('minimize 1e-10 1e-10 1000000 10000000\n')
        o.write('run 1\n')
        o.write('write_data res.dat\n')
    target = subprocess.check_output([LAMMPS,'-in','script.txt'],universal_newlines=True).split('\n')
    o=enumerate(target)
    for idx, line in o:
        tmp = line.split()
        if len(tmp) > 0:
            if tmp[0] =='Step':
                next(o)
                step, T, PE = list(map(float,next(o)[-1].split()[:3]))
    return PE

def set_pot(FILE,pot=pot,pot_type=pot_type):
    if pot_type == 'tersoff':
        FILE.write('pair_style  tersoff\n')
        FILE.write(f'pair_coeff  * * "{pot}" Si \n')
    elif pot_type == 'tersoff/mod':
        FILE.write('pair_style  tersoff/mod\n')
        FILE.write(f'pair_coeff  * * "{pot}" Si \n')
    elif pot_type == 'sw':
        FILE.write('pair_style  sw\n')
        FILE.write(f'pair_coeff  * * "{pot}" Si \n')
    return 0

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
