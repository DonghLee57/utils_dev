import sys
import numpy as np
import matplotlib.pyplot as plt
from ase.formula import Formula
import ase.io

def main():
    plot_pdos = 0
    NORBIT = 9
    struct = ase.io.read('POSCAR',format='vasp')
    symbols = struct.get_chemical_symbols()
    ntypes = Formula(str(struct.symbols)).count()
    
    MYDOS = DOS(NORBIT)
    MYDOS.read(sys.argv[1])
    
    ENERGY = np.linspace(MYDOS.EMIN, MYDOS.EMAX, MYDOS.NEDOS)
    PDOS = {}
    ATOMS = {}
    for e in ntypes:
        ATOMS[e] = [] 
        for i, item in enumerate(symbols):
            if e==item: ATOMS[e].append(i)


    fig, ax = plt.subplots()
    c_index = list(ntypes.keys())
    ax.plot(ENERGY-MYDOS.fermi, MYDOS.TDOS[:,1],'k')
    if MYDOS.ISPIN == 2:
        ax.plot(ENERGY-MYDOS.fermi, -MYDOS.TDOS[:,2],'k')
    else:
        ax.set_ylim(bottom=0)

    if plot_pdos:
        if MYDOS.ISPIN == 2:
            ORBIT  = np.arange(1,NORBIT*MYDOS.ISPIN+1,2)
            ORBIT2 = np.arange(2,NORBIT*MYDOS.ISPIN+1,2)
            for e in ntypes:
                PDOS[e] = np.sum(np.sum(MYDOS.DOS[ATOMS[e]].T[ORBIT],axis=0),axis=1)
                ax.plot(ENERGY-MYDOS.fermi, PDOS[e])
            for e in ntypes:
                PDOS[e] = np.sum(np.sum(MYDOS.DOS[ATOMS[e]].T[ORBIT2],axis=0),axis=1)
                ax.plot(ENERGY-MYDOS.fermi, -PDOS[e],c='C'+str(c_index.index(e)))

        else:
            ORBIT  = np.arange(1,NORBIT*MYDOS.ISPIN+1) 
            for e in ntypes:
                PDOS[e] = np.zeros(MYDOS.NEDOS) 
                for i in range(len(ATOMS[e])):
                    PDOS[e] = np.sum(MYDOS.DOS[ATOMS[e]].T[ORBIT,],axis=0)
                ax.plot(ENERGY-MYDOS.fermi, PDOS[e],c='C'+str(c_index.index(e)))
    ax.axvline(0,c='gray', ls='--')
    ax.set_xlim([-2, 2])
    ax.set_xlabel(r'E-E$_f$ (eV)',fontsize=12)
    ax.set_ylabel(r'DOS (a.u.)',fontsize=12)
    ax.set_yticks([])
    plt.tight_layout()
    plt.show()

#--------------------------------------------------------------------------------------
class DOS:
    def __init__(self,NORBIT=9):
        self.NIONS = 0
        self.EMAX, self.EMIN, self.NEDOS, self.fermi = 0, 0, 0, 0
        self.ISPIN = 1
        self.NORBIT= NORBIT
        self.DOS = []
        self.TDOS= []

    def read(self, DOSCAR):
        with open(DOSCAR,'r') as o: tmp = enumerate(o.readlines())
        self.NIONS = int(next(tmp)[1].split()[0])
        for i in range(4): next(tmp)
        self.EMAX, self.EMIN, self.NEDOS, self.fermi, _ =  map(float,next(tmp)[1].split())
        self.NEDOS = int(self.NEDOS)

        for e in range(self.NEDOS):
            _, dat = next(tmp)
            if len(self.TDOS) == 0:
                if len(dat.split()) <=3: self.ISPIN = 1
                else: self.ISPIN = 2 
                self.TDOS = np.zeros((self.NEDOS,self.ISPIN*2+1))
                self.DOS  = np.zeros((self.NIONS, self.NEDOS, self.ISPIN*self.NORBIT+1))
            else:
                self.TDOS[e] += np.array(list(map(float,dat.split())))
        n = 0
        for idx, line in tmp:
            for e in range(self.NEDOS):
                _, dat = next(tmp)
                self.DOS[n][e] += np.array(list(map(float,dat.split())))
            n+=1
        return 0

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
