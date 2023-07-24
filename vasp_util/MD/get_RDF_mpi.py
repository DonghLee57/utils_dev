import sys
import numpy as np
from ase import Atoms
from ase.io import read
import matplotlib.pyplot as plt
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

r_max= 6
dr = 0.01
tag = ['rdf','prdf']

start = int(sys.argv[2])
end = int(sys.argv[3])
step = int(sys.argv[4])

if len(sys.argv[:]) < 3:
    base = read(sys.argv[1],index=':')
elif len(sys.argv[:]) ==3:
    base = read(sys.argv[1],index=start)
elif len(sys.argv[:]) > 3:
    base = read(sys.argv[1],index=slice(start,end,step))
nimg = len(base)
volume = np.fabs(np.linalg.det(base[0].cell))
symbols = np.array(base[0].get_chemical_symbols())
types = list(set(symbols))

if 'rdf' in tag: 
    bins = np.arange(dr/2, r_max, dr)
    RDF = np.zeros((nimg,len(bins)-1))
    for idx in range(nimg):
        positions = base[idx].positions
        NIONS = len(positions)
        dist = np.zeros((NIONS,NIONS))
        psize = NIONS//size
        if rank == size-1:
            for a in range(rank*psize, NIONS):
                dist[a] = base[idx].get_distances(a,range(NIONS),mic=True)
        else:
            for a in range(rank*psize,(rank+1)*psize):
                dist[a] = base[idx].get_distances(a,range(NIONS),mic=True)
        dist = comm.reduce(dist, MPI.SUM)
        dist = comm.bcast(dist, root=0)
        res = np.histogram(dist,bins=bins)
        rdf = volume*res[0]/(NIONS**2*4*np.pi*res[1][:-1]**2*dr)
        RDF[idx] = rdf
    if rank==0:
        RDF = np.mean(np.array(RDF),axis=0)
        fig,ax = plt.subplots()
        ax.plot(res[1][:-1], RDF)
        print(f"Peak height: {max(RDF):.2f}")
        print(f"Peak position: {res[1][np.argmax(RDF)]:.2f}")
        ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=12)
        ax.set_ylabel('g(r)',fontsize=12)
        ax.set_xlim([0,r_max])
        ax.set_ylim(bottom=0)
        plt.savefig('Total_RDF.png')
        plt.clf()

if 'prdf' in tag:
    if rank == 0:
        combinations = []
        for idx, i in enumerate(types):
            for jdx, j in enumerate(types):
                if idx <= jdx:
                    combinations.append((i,j))
        
        bins = np.arange(dr/2, r_max, dr)
        fig,ax = plt.subplots()
    else:
        combinations = None
        bins = np.array([])
    combinations = comm.bcast(combinations, root=0)
    bins = comm.bcast(bins, root=0)
   
    for pidx, pair in enumerate(combinations):
        i = np.where(symbols == pair[0])[0]
        j = np.where(symbols == pair[1])[0]
        RDF = np.zeros((nimg,len(bins)-1))
        for idx in range(nimg):
            positions = base[idx].positions
            dist = np.zeros((len(i),len(j)))
            psize = len(i)//size
            if rank == size-1:
                for a in range(rank*psize,len(i)):
                    dist[a] = base[idx].get_distances(i[a],j,mic=True)
            else:
                for a in range(rank*psize,(rank+1)*psize):
                    dist[a] = base[idx].get_distances(i[a],j,mic=True)
            dist = comm.reduce(dist,op=MPI.SUM,root=0)
            dist = comm.bcast(dist, root=0)
            res = np.histogram(dist,bins=bins)
            rdf = volume*res[0]/(len(i)*len(j)*4*np.pi*res[1][:-1]**2*dr)
            RDF[idx] = rdf
        if rank == 0:
            RDF = np.mean(np.array(RDF),axis=0)
            ax.plot(res[1][:-1], RDF,label=f'{pair[0]}-{pair[1]}')
            print(f"{pair[0]}-{pair[1]}")
            print(f"Peak height: {max(RDF):.2f}")
            print(f"Peak position: {res[1][np.argmax(RDF)]:.2f}")
    if rank==0:
        ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=12)
        ax.set_ylabel('g(r)',fontsize=12)
        ax.set_xlim([0,r_max])
        ax.set_ylim(bottom=0)
        plt.legend()
        plt.savefig('Partial_RDF.png')
        plt.clf()
