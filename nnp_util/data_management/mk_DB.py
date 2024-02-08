import sys
import os, glob
import sqlite3

PRE_PATH =  '/tmp/'
DB_SAVE  = f'{PRE_PATH}/Lee/my.db'
MEMBERS  = ['Lee']
MAXREADOUTCAR = 1000

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

# Main - Construct Database
def main():
    global DB_SAVE, MAT_MEMBERS, PRE_PATH
    sql = sqlite3.connect(DB_SAVE)
    cur = sql.cursor()
    cur.execute('''
    CREATE TABLE mydb (
    Elements TEXT,
    Structure_type TEXT,
    File_path TEXT,
    init TEXT,
    end  TEXT,
    step TEXT,
    comment TEXT)''')
    for user in MEMBERS:
        path_user_db = f'{PRE_PATH}/{user}/DB/'
        strlist = glob.glob(f'{path_user_db}/*')
        for FILE in strlist:
            print(FILE)
            parse_insert(cur, FILE)
    sql.commit()
    sql.close()
    return 0

# Functions
def parse_insert(CURSOR, FILE):
    # parsing data in structure_list format from SIMPLE-NN
    global EDICT
    with open(FILE,'r') as o: tmp = o.readlines()
    comment = 'No comment'
    for idx in range(len(tmp)):
        line = tmp[idx]
        line = line.strip()
        if len(line) > 0 :
                if line[0] == '#':
                    comment = ' '.join(line.split()[1:])
                elif line[0] == '[':
                    structure_type = ''.join(line.split()).strip('[]')
                else:
                    if len(line.split()) == 2 and len(line.split(':')) == 3:
                        [path, slice_idx] = line.split()
                        # : split X case ?? 
                        [init, end, step] = slice_idx.split(':')
                        if os.path.isfile(path):
                            elements = collect_elem_from_outcar(path)     
                            CURSOR.execute('INSERT INTO mydb (Elements, Structure_type, File_path, init, end, step, comment) VALUES (?, ?, ?, ?, ?, ?, ?)',
                                          (elements, structure_type, path, init, end, step, comment))
                        else:
                            print(f"No OUTCAR: {path}")
                    else:
                        print(f"Illegal format: {line}" )

def collect_elem_from_outcar(OUTCAR):
    global MAXREADOUTCAR
    elem = []
    with open(OUTCAR, 'r') as o: tmp = o.readlines()[:MAXREADOUTCAR]
    for i in range(len(tmp)):
        line = tmp[i].split()
        if 'POTCAR:' in line:
            elem.append(line[-2].split('_')[0])
    elem = list(set(elem))
    elem.sort(key=lambda x: EDICT[x])
    elements = '-'.join(elem)
    return elements

#
if __name__ == "__main__":
    main()
