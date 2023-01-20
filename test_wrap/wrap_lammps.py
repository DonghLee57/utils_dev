myLAMMPS = '/path/lammps/python'
import sys
sys.path.append(myLAMMPS)

import lammps
import ase.io
from ase import Atoms

def main():
    # For calculating energy difference, energies of the reactant, product molecules should calculated.
    mol_list = glob.glob('./molecules/*')
    mol_E = {}
    for idx, m in enumerate(mol_list):
        mol = ase.io.read(m,index=':',format='vasp')
        types = list(set(mol[0].get_chemical_symbols()))
        ase.io.write(str4lmp, images=mol, format='lammps-data')
        write_lmp_script(str4lmp,nnp,types)

        lmp = lammps.lammps()
        lmp.file('py_script.in')
        mol_E[m.split('_')[-1]] = get_lmp_log('energy')
        lmp.close()

    input_str = sys.argv[1]
    str4lmp = sys.argv[2]
    nnp = sys.argv[3]
    
    info = readPOSCAR(input_str)
    poscar = ase.io.read(input_str,index=':',format='vasp')

    ase.io.write(str4lmp, images=poscar, format='lammps-data')
    write_lmp_script(str4lmp,nnp,info[1])

    lmp = lammps.lammps()
    lmp.file('py_script.in')
    E0 = get_lmp_log('energy')
    lmp.close()

def write_lmp_script(str_name, pot_name, elements):
    periodic_table = [
    'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn','Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd','In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba',    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
               'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi','Po', 'At', 'Rn',
    'Fr', 'Ra',    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk','Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
               'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc','Lv', 'Ts', 'Og']
    e_dict ={}
    for Z,symbol in enumerate(periodic_table):
        e_dict[symbol] = Z+1
    # In the present ASE, num(Cl) < num(N)...
    e_dict['Cl'] = 6
    mass_dict={'H':1, 'N':14.001, 'O':15.999, 'Cl':35.453,'Ti':47.880}
    o=open(pot_name,'r').readline()
    pot_elements = o.split()[1:]
    o=open('py_script.in','w')
    header='units           metal\n'
    header+='boundary        p p p\n'
    header+='atom_style      atomic\n'
    header+='box tilt        large\n'
    header+='read_data       %s\n' %str_name
    header+='\npair_style      nn\n'

    elements.sort(key=lambda x: e_dict[x])
    for i in range(len(elements)):
        pot_elements.remove(elements[i])
    coeff_ele = elements + pot_elements

    header+='pair_coeff      * * %s ' %pot_name
    header+=' '.join(coeff_ele)
    header+='\n\n'
    for idx, e in enumerate(coeff_ele):
        header+=f'mass\t{idx+1} {mass_dict[e]}\n'

    header+='\nrun 0'
    o.write(header)
    return 1


def get_lmp_log(prop):
    o=enumerate(open('log.lammps'))
    for idx, line in o:
        tmp = line.split()
        if len(tmp) > 0:
            if tmp[0] =='Step':
                step, T, PE = list(map(float,next(o)[-1].split()[:3]))
    if prop.lower() == 'energy': return PE

##########################
if __name__ == "__main__":
    main()
