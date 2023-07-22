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

base = read(sys.argv[1],format='vasp')
volume = np.fabs(np.linalg.det(base.cell))
positions = base.positions.copy()
symbols = np.array(base.get_chemical_symbols())
types = list(set(symbols))

if 'rdf' in tag: 
    bins = np.arange(dr/2, r_max, dr)
    dist = []
    NIONS = len(positions)
    psize = NIONS//size
    if rank == size-1:
        for a in range(rank*psize,NIONS):
            for b in range(NIONS):
                if a != b:
                    dist.append(base.get_distance(a,b,mic=True))
    else:
        for a in range(rank*psize,(rank+1)*psize):
            for b in range(NIONS):
                if a != b:
                    dist.append(base.get_distance(a,b,mic=True))
    dist = comm.gather(dist, root=0)
    if rank == 0:
        res = np.histogram(dist,bins=bins)
        fig,ax = plt.subplots()
        ax.plot(res[1][:-1],volume*res[0]/(NIONS**2*4*np.pi*res[1][:-1]**2*dr))
        ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=12)
        ax.set_ylabel('g(r)',fontsize=12)
        ax.set_xlim([0,r_max])
        ax.set_ylim(bottom=0)
        plt.savefig('Total_RDF.png')
        plt.clf()

if 'prdf' in tag: 
    combinations = []
    for idx, i in enumerate(types):
        for jdx, j in enumerate(types):
            if idx <= jdx:
                combinations.append((i,j))
    
    bins = np.arange(dr/2, r_max, dr)
    fig,ax = plt.subplots()
    for idx, pair in enumerate(combinations):
        i = np.where(symbols == pair[0])[0]
        j = np.where(symbols == pair[1])[0]
        dist = []
        psize = len(i)//size
        if rank == size-1:
            for a in range(rank*psize,len(i)):
                for b in j:
                    if a != b:
                        dist.append(base.get_distance(i[a],b,mic=True))
        else:
            for a in range(rank*psize,(rank+1)*psize):
                for b in j:
                    if a != b:
                        dist.append(base.get_distance(i[a],b,mic=True))
        dist = comm.gather(dist, root=0)
        if rank == 0:
            res = np.histogram(dist,bins=bins)
            plt.plot(res[1][:-1],volume*res[0]/(len(i)*len(j)*4*np.pi*res[1][:-1]**2*dr),label=f'{pair[0]}-{pair[1]}')
    if rank == 0:
        ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=12)
        ax.set_ylabel('g(r)',fontsize=12)
        ax.set_xlim([0,r_max])
        ax.set_ylim(bottom=0)
        plt.legend()
        plt.savefig('Partial_RDF.png')
        plt.clf()
