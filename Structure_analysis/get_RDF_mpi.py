import sys
import numpy as np
from ase import Atoms
from ase.io import read
import matplotlib.pyplot as plt
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

### Pre-defined parameters
r_max= 6
dr = 0.01

###
def main():
    # For a single image
    if True:
        tmp = read("POSCAR", format="vasp")
        symbols = np.array(tmp.get_chemical_symbols())
        types = list(set(symbols))
        RDF = get_rdf(tmp)
        PRDF = get_prdf(tmp, types)
        
    # For temporal averaging
    if False:
        tmp = read("XDATCAR", index=slice(start,end,step), format="vasp-xdatcar")
        nimg = len(tmp)
        RDF = None
        #PRDF = None
        for idx, img in enumerate(tmp):
            symbols = np.array(img.get_chemical_symbols())
            types = list(set(symbols))
            if True:
                if RDF == None: RDF = get_rdf(img)
                else:           RDF[1] += get_rdf(img)[1]            
            if True:
                if PRDF == None: PRDF = get_prdf(img, types)
                else:
                    comb = PRDF[1].keys()
                    for idx, pair in enumerate(comb):
                        PRDF[1][pair] += get_prdf(img, types)[1][pair]
        for idx, pair in enumerate(comb): PRDF[1][pair] /= nimg
        
    # Plotting data
    if False:
        plot_rdf(RDF)
        plot_prdf(PRDF)

###
def get_rdf(Obj):
    bins = np.arange(dr/2, r_max, dr)
    positions = Obj.positions
    volume = np.fabs(np.linalg.det(Obj.cell))
    NIONS = len(positions)
    dist = np.zeros((NIONS,NIONS))
    psize = NIONS//size
    if rank == size-1:
        for a in range(rank*psize, NIONS):
            dist[a] = Obj.get_distances(a,range(NIONS),mic=True)
    else:
        for a in range(rank*psize,(rank+1)*psize):
            dist[a] = Obj.get_distances(a,range(NIONS),mic=True)
    dist = comm.reduce(dist, MPI.SUM)
    dist = comm.bcast(dist, root=0)
    res = np.histogram(dist,bins=bins)
    rdf = volume*res[0]/(NIONS**2*4*np.pi*res[1][:-1]**2*dr)
    return res[1][:-1], rdf

def plot_rdf(RDF):
    fig,ax = plt.subplots()
    ax.plot(RDF[0], RDF[1])
    print(f"Peak height: {np.max(RDF[1]):.2f}")
    print(f"Peak position: {RDF[0][np.argmax(RDF[1])]:.2f}")
    ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=12)
    ax.set_ylabel('g(r)',fontsize=12)
    ax.set_xlim([0,r_max])
    ax.set_ylim(bottom=0)
    plt.savefig('Total_RDF.png')
    return 0

def get_prdf(Obj, types, targets=None):
    # targets: ("element1", "element2")
    # ex) ("Si", "O")
    bins = np.arange(dr/2, r_max, dr)
    volume = np.fabs(np.linalg.det(Obj.cell))
    if targets == None:
        if rank == 0:
            combinations = []
            for idx, i in enumerate(types):
                for jdx, j in enumerate(types):
                    if idx <= jdx:
                        combinations.append((i,j))
            fig,ax = plt.subplots()
        else:
            combinations = None
        combinations = comm.bcast(combinations, root=0)
    else:
        combinations = [targets]       
    prdf = {}
    for pidx, pair in enumerate(combinations):
        i = np.where(symbols == pair[0])[0]
        j = np.where(symbols == pair[1])[0]
        for idx in range(nimg):
            positions = Obj.positions
            dist = np.zeros((len(i),len(j)))
            psize = len(i)//size
            if rank == size-1:
                for a in range(rank*psize,len(i)):
                    dist[a] = Obj.get_distances(i[a],j,mic=True)
            else:
                for a in range(rank*psize,(rank+1)*psize):
                    dist[a] = Obj.get_distances(i[a],j,mic=True)
            dist = comm.reduce(dist,op=MPI.SUM,root=0)
            dist = comm.bcast(dist, root=0)
            res = np.histogram(dist,bins=bins)
            rdf = volume*res[0]/(len(i)*len(j)*4*np.pi*res[1][:-1]**2*dr)
            prdf[pair] = rdf.copy()
    return res[1][:-1], prdf

def plot_prdf(PRDF):
    fig,ax = plt.subplots()
    comb = PRDF[1].keys()
    for idx, pair in enumerate(comb):
        ax.plot(PRDF[0], PRDF[1][pair], label=f'{pair}')
        print(f"{pair}")
        print(f"Peak height: {np.max(PRDF[1][pair]):.2f}")
        print(f"Peak position: {PRDF[0][np.argmax(PRDF[1][pair])]:.2f}")
    ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=12)
    ax.set_ylabel('g(r)',fontsize=12)
    ax.set_xlim([0,r_max])
    ax.set_ylim(bottom=0)
    ax.legend()
    plt.savefig('Partial_RDF.png')
    return 0

###
if __name__ == "__main__":
    main()
