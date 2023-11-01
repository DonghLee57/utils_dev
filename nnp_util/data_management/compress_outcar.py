import sys
#def compress_outcar(filename):
"""
    Compress VASP OUTCAR file for fast file-reading in ASE.
    Compressed file (tmp_comp_OUTCAR) is temporarily created in the current directory.

    :param str filename: filename of OUTCAR

    supported properties:

    - atom types
    - lattice vector(cell)
    - free energy
    - force
    - stress
"""
FILE = sys.argv[1]
OUT  = sys.argv[2]

with open(FILE, 'r') as fil, open(OUT, 'w') as res:
    minus_tag = 0
    line_tag = 0
    ions_key = 0
    for line in fil:
        if 'POTCAR:' in line:
            res.write(line)
        if 'POSCAR:' in line:
            res.write(line)
        elif 'ions per type' in line and ions_key == 0:
            res.write(line)
            ions_key = 1
        elif 'direct lattice vectors' in line:
            res.write(line)
            minus_tag = 3
        elif 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM' in line:
            res.write(line)
            minus_tag = 4
        elif 'POSITION          ' in line:
            res.write(line)
            line_tag = 3
        elif 'FORCE on cell =-STRESS' in line:
            res.write(line)
            minus_tag = 15
        elif 'Iteration' in line:
            res.write(line)
        elif minus_tag > 0:
            res.write(line)
            minus_tag -= 1
        elif line_tag > 0:
            res.write(line)
            if '-------------------' in line:
                line_tag -= 1
