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
# RDF-related parameter
r_max= 6
dr = 0.01
# ADF-related parameter
angle_lim = [0, 180]
dangle=1.0
# plotting style
fs=12

###
def main():
    # For a single image
    if True:
        tmp = read("POSCAR", format="vasp")
        symbols = np.array(tmp.get_chemical_symbols())
        types = list(set(symbols))
        RDF = get_rdf(tmp)
        PRDF = get_prdf(tmp, types)
        ADF = get_adf(tmp, ["B","A","A"], [r1, r2], angle_lim=angle_lim):
        
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
        plot_adf(ADF, triplet_label='A-B-A')

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
    ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=fs)
    ax.set_ylabel('g(r)',fontsize=fs)
    ax.set_xlim([0,r_max])
    ax.set_ylim(bottom=0)
    plt.savefig('Total_RDF.png')
    return 0

def get_prdf(Obj, types, targets=None):
    # targets: ["element1", "element2"]
    # ex) ["Si", "O"]
    bins = np.arange(dr/2, r_max, dr)
    volume = np.fabs(np.linalg.det(Obj.cell))
    if targets == None:
        if rank == 0:
            combinations = []
            for idx, i in enumerate(types):
                for jdx, j in enumerate(types):
                    if idx <= jdx:
                        combinations.append([i,j])
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
    ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=fs)
    ax.set_ylabel('g(r)',fontsize=fs)
    ax.set_xlim([0,r_max])
    ax.set_ylim(bottom=0)
    ax.legend()
    plt.savefig('Partial_RDF.png')
    return 0
    
def get_adf(Obj, targets, cutoff, angle_lim=[0, 180], expr='degree'):
    # To facilitate a clear bond angle analysis, the target triplet with cutoff distances should be provided.
    # targets: ["Center", "Neighbor_1", "Neighbor_2"]
    # ex) ["Si", "Si", "O"]
    # cutoff: [r_1, r_2]
    # ex) [2.6, 2.0]
    bins = np.arange(angle_lim[0], angle_lim[1]+0.001, dangle)
    if expr != 'degree': bins = bins*np.pi/180
    symbols = np.array(tmp.get_chemical_symbols())
    cidx = np.where( targets[0] == symbols )[0]
    nidx = np.where( targets[1] == symbols )[0]
    midx = np.where( targets[2] == symbols )[0]
    theta = []
    psize = len(cidx)//size
    if rank == size-1:
        for c in range(rank*psize,len(cidx)):
            for n in range(len(nidx):
                vec1 = Obj.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                dist1 = np.linalg.norm(vec)
                if dist1 < cutoff[0]:
                    for m in range(len(midx)):
                        vec2 = Obj.get_distances(cidx[c], midx[m], mic=True, vector=True)
                        dist2 = np.linalg.norm(vec)
                        if dist2 < cutoff[1]:
                            theta.append(np.arccos(np.dot(vec1, vec2)/dist1/dist2))
    else:
        for c in range(rank*psize,(rank+1)*psize):
            for n in range(len(nidx):
                vec1 = Obj.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                dist1 = np.linalg.norm(vec)
                if dist1 < cutoff[0]:
                    for m in range(len(midx)):
                        vec2 = Obj.get_distances(cidx[c], midx[m], mic=True, vector=True)
                        dist2 = np.linalg.norm(vec)
                        if dist2 < cutoff[1]:
                            theta.append(np.arccos(np.dot(vec1, vec2)/dist1/dist2))
    theta = comm.reduce(theta,op=MPI.SUM,root=0)
    theta = comm.bcast(theta, root=0)
    res = np.histogram(theta,bins=bins)
    if expr != 'degree': return res[1][:-1], res[0]
    else:                return res[1][:-1], res[0]/np.pi*180
    
 def plot_adf(ADF, expr='degree', triplet_label=''):
    fig,ax = plt.subplots()
    if triplet_label != '': fig.suptitle(triplet_label, fontsize=fs*1.5)
    ax.plot(ADF[0], ADF[1])
    if expr != 'degree': ax.set_xlabel(r'Bond angle (rad)',fontsize=fs)
    else:                ax.set_xlabel(r'Bond angle ($^\circ$)',fontsize=fs)
    ax.set_ylabel('Counts', fontsize=fs)
    ax.set_xlim([ADF[0][0], ADF[0][-1]])
    ax.set_ylim(bottom=0)
    plt.savefig(f'ADF_{triplet_label}.png')
    return 0 

###
if __name__ == "__main__":
    main()
