# python3.X
import sys

#--------------------------------------------------------------------------------------
def main():
    FILE = sys.argv[1]
    FLAG_It, FLAG_scf = check_e_conv(FILE)
    print(f'{"FILENAME"} {"OUTCAR_Writing"} {"Iteration"} {"SCF"}')
    if not check_outcar(FILE) or not FLAG_scf or not FLAG_It:
        print(f'{FILE} {check_outcar(FILE)} {FLAG_It} {FLAG_scf}')

#--------------------------------------------------------------------------------------
def check_outcar(PATH):
    FLAG_OUT = False
    # print(path)
    with open(PATH,'r') as o: tmp = o.readlines()[-1]
    if len(tmp.split()) > 0:
        if 'Voluntary' == tmp.split()[0]: FLAG_OUT = True
        else: FLAG_OUT = False
    return FLAG_OUT

def check_e_conv(PATH):
    FLAG_Econv = True
    FLAG_It = False
    with open(PATH,'r') as o: tmp = o.readlines()
    for i in range(len(tmp)):
        line = tmp[i].split()
        if 'NELM' in line: NELM = line[2][:-1]
        if 'NSW'  in line:  NSW = line[2]
        if 'Iteration' in line:
            if NELM == line[3][:-1]:
                FLAG_Econv = False
            if NSW == line[2][:-1]:
                FLAG_It = True

    return FLAG_It,FLAG_Econv
#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
