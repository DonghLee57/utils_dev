# python3.X
import sys
import os
import glob

#--------------------------------------------------------------------------------------
def main():
    check_outcar(sys.argv[1])

#--------------------------------------------------------------------------------------
def check_outcar(PATH):
    if not os.path.isfile(PATH):
        print(f"!!! NO SUCH FILE: {PATH} !!!")
        return 0
    else:
        with open(PATH,'r') as o: tmp = o.readlines()
        # Normal END?
        if 'Voluntary' != tmp[-1].split()[0]:
            print(f"!!! ABNORMAL END: {PATH} !!!")
            return 0
        # Electronic step Convergence ?
        for i in range(len(tmp)):
            line = tmp[i].split()
            if 'NELM' in line: NELM = line[2][:-1]
            if 'Iteration' in line:
                if NELM == line[3][:-1]:
                    print(f"!!! SCF NOT CONVERGED: {PATH} !!!")
                    return 0
#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
