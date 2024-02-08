# python3.x
# ver.1 @24.02.06
#
import sys
import os
import sqlite3
import itertools

PRE_PATH =  '/tmp/'
DB_SAVE  = f'{PRE_PATH}/Lee/my.db'

#from ase.data import atomic_numbers
EDICT = {'X': 0, 'H': 1, 'He': 2,\
         'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,\
         'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,\
         'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,\
         'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,\
         'Cs': 55, 'Ba': 56,\
         'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,\
         'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,\
         'Fr': 87, 'Ra': 88,\
         'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,\
         'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

# Main
def main():
    global DB_SAVE, EDICT
    elements = sys.argv[1:]
    elements = list(set(elements))
    elements.sort(key=lambda x: EDICT[x])

    sql = sqlite3.connect(DB_SAVE)
    cur = sql.cursor()
    for N in range(len(elements)):
        combination = list(itertools.combinations(elements,N+1))
        for my_combi in combination:
            # Search database
            my_combi = list(my_combi)
            find = "-".join(my_combi)
            cur.execute(f'SELECT DISTINCT Elements, Structure_type, File_path, init, end, step FROM MLPDB WHERE Elements LIKE "{find}" ORDER BY Elements ASC, Structure_type ASC')
            res = cur.fetchall()

            # Generate structure list
            if len(res) > 0:
                write_structure_list(res)
            else:
                print(f"Data consisting of {find} not in my.db")
    sql.commit()
    sql.close()
    return 0

def write_structure_list(info):
    o = open(f'./DB_{info[0][0]}','w')
    memo_stype = None
    for tid, line in enumerate(info):
        stype = line[:2]
        path, init, end, step = line[2:]
        if memo_stype != stype:
            memo_stype = stype
            o.write(f'\n[{"_".join(stype)} : 1]\n')
        o.write(f'{path} {init}:{end}:{step}\n')
    o.close()

#
if __name__ == "__main__":
    main()
