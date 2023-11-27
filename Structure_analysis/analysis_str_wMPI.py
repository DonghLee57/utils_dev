import sys, os
import numpy as np
from ase import Atoms
from ase.io import read
import matplotlib.pyplot as plt
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

### Pre-defined parameters
start, end, step = 0, 50 , 1
# RDF-related parameter
r_max= 3
dr = 0.01
RDF_FILE = 'rdf.out'
PRDF_FILE = 'prdf.out'
# ADF-related parameter
triplet    =  ["Si","O","O"]  #["Center", "Neighbor_1", "Neighbor_2"]
rc4angle   =  [1.6 , 1.6]     # Pair of bond lengths for calculating angle
angle_lim  =  [0, 180]        # Minimum and maximum of bin
dangle     =  1.0             # bin size in histgram
# Diffusivity
SELECT_ATOMS = np.arange(0,2)
MSD_FILE = 'msd.out'

# plotting style
fs=12

###
def main():
    # For a single image
    if False:
        tmp = read(sys.argv[1], format="vasp")
        symbols = np.array(tmp.get_chemical_symbols())
        types = list(set(symbols))
        #RDF = get_rdf(tmp, savefile=RDF_FILE)
        #PRDF = get_prdf(tmp, types, savefile=PRDF_FILE)
        ADF = get_adf(tmp, triplet, rc4angle, angle_lim=angle_lim)
        
    # For temporal averaging
    if True:
        tmp = read("XDATCAR", index=slice(start,end,step), format="vasp-xdatcar")
        #tmp = read(sys.argv[1], index=slice(start,end,step), format="lammps-dump-text")
        nimg = len(tmp)
        if False:
            RDF = None
            for idx, img in enumerate(tmp):
                symbols = np.array(img.get_chemical_symbols())
                types = list(set(symbols))
                if RDF == None: RDF     = get_rdf(img)
                else:           RDF[1] += get_rdf(img)[1]
            RDF[1] /= nimg
        if False:
            PRDF = None
            for idx, img in enumerate(tmp):
                symbols = np.array(img.get_chemical_symbols())
                types = list(set(symbols))
                if PRDF == None: PRDF = get_prdf(img, types)
                else:
                    comb = PRDF[1].keys()
                    for idx, pair in enumerate(comb):
                        PRDF[1][pair] += get_prdf(img, types)[1][pair]
            for idx, pair in enumerate(comb):
                PRDF[1][pair] /= nimg
        if False:
            ADF = None
            for idx, img in enumerate(tmp):
                symbols = np.array(img.get_chemical_symbols())
                types = list(set(symbols))
                if ADF == None:  ADF     = get_adf(img, triplet, rc4angle, angle_lim=angle_lim)
                else:            ADF[1] += get_adf(img, triplet, rc4angle, angle_lim=angle_lim)
            ADF /= nimg
        if True:
            UNWRAP = unwrapping(tmp, targets=SELECT_ATOMS)
            # One-line
            MSD = np.sum(np.linalg.norm(UNWRAP - UNWRAP[0],axis=-1),axis=-1)/len(SELECT_ATOMS)
            # Parallel?
            #MSD = ??

    # Plotting data
    if True and rank == 0:
        #if os.path.isfile(RDF_FILE): RDF = np.loadtxt(RDF_FILE)
        #plot_rdf(RDF)
        #plot_prdf(PRDF)
        #plot_adf(ADF, triplet_label=f'(triplet[1])-{triplet[0]}-{triplet[2]}')
        plot_msd(MSD, unit='ps')
        pass

###
def get_rdf(Obj, savefile='rdf.out'):
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
    if rank == 0: np.savetxt(savefile, [res[1][:-1], rdf], fmt='%.4f')
    return res[1][:-1], rdf

def plot_rdf(RDF, savefile='Total_RDF.png'):
    fig,ax = plt.subplots()
    ax.plot(RDF[0], RDF[1])
    print(f"Peak height: {np.max(RDF[1]):.2f}")
    print(f"Peak position: {RDF[0][np.argmax(RDF[1])]:.2f}")
    ax.set_xlabel(r'Distance ($\mathrm{\AA}$)',fontsize=fs)
    ax.set_ylabel('g(r)',fontsize=fs)
    ax.set_xlim([0,r_max])
    ax.set_ylim(bottom=0)
    plt.savefig(savefile)
    return 0

def get_prdf(Obj, types, targets=None, savefile='prdf.out'):
    # targets: ["element1", "element2"]
    # ex) ["Si", "O"]
    bins = np.arange(dr/2, r_max, dr)
    symbols = np.array(Obj.get_chemical_symbols())
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
        key = f'{pair[0]}-{pair[1]}'
        prdf[key] = rdf.copy()
        if rank == 0:
            np.savetxt(key+'_'+savefile, [res[1][:-1], prdf[key]], fmt='%.4f')
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
    symbols = np.array(Obj.get_chemical_symbols())
    cidx = np.where( targets[0] == symbols )[0]
    nidx = np.where( targets[1] == symbols )[0]
    midx = np.where( targets[2] == symbols )[0]
    theta = []
    psize = len(cidx)//size
    if targets[1] == targets[2]:
        if rank == size-1:
            for c in range(rank*psize,len(cidx)):
                for n in range(len(nidx)):
                    vec1 = Obj.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                    dist1 = np.linalg.norm(vec1)
                    if dist1 < cutoff[0]:
                        for m in range(len(midx)):
                            if m > n:
                                vec2 = Obj.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                dist2 = np.linalg.norm(vec2)
                                if dist2 < cutoff[1]:
                                    theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
        else:
            for c in range(rank*psize,(rank+1)*psize):
                for n in range(len(nidx)):
                    vec1 = Obj.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                    dist1 = np.linalg.norm(vec1)
                    if dist1 < cutoff[0]:
                        for m in range(len(midx)):
                            if m > n:
                                vec2 = Obj.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                dist2 = np.linalg.norm(vec2)
                                if dist2 < cutoff[1]:
                                    theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
    else:
        if rank == size-1:
            for c in range(rank*psize,len(cidx)):
                for n in range(len(nidx)):
                    vec1 = Obj.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                    dist1 = np.linalg.norm(vec1)
                    if dist1 < cutoff[0]:
                        for m in range(len(midx)):
                            if m != n:
                                vec2 = Obj.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                dist2 = np.linalg.norm(vec2)
                                if dist2 < cutoff[1]:
                                    theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
        else:
            for c in range(rank*psize,(rank+1)*psize):
                for n in range(len(nidx)):
                    vec1 = Obj.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                    dist1 = np.linalg.norm(vec1)
                    if dist1 < cutoff[0]:
                        for m in range(len(midx)):
                            if m != n:
                                vec2 = Obj.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                dist2 = np.linalg.norm(vec2)
                                if dist2 < cutoff[1]:
                                    theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
    theta = comm.reduce(theta,op=MPI.SUM,root=0)
    theta = comm.bcast(theta, root=0)
    if expr == 'degree': theta = np.array(theta)*180/np.pi
    res = np.histogram(theta,bins=bins)
    return res[1][:-1], res[0]
    
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

def unwrapping(Obj, targets):
    IMG0 = Obj[0].get_positions()[targets]
    LAT = Obj[0].cell.copy()
    OLD = Obj[0].get_scaled_positions()
    UNWRAP_TRJ = np.zeros((len(Obj), len(targets), 3))
    psize = len(targets)//size
    if rank == size-1:
        parts = np.arange(rank*psize,len(targets))
        for idx in range(1, len(Obj)):
            NEW = Obj[idx].get_scaled_positions()
            diff = NEW[parts] - OLD[parts]
            UNWRAP_TRJ[idx][parts] = UNWRAP_TRJ[idx-1][parts] + (np.round(-diff)+diff)@LAT 
            OLD = Obj[idx].get_scaled_positions()
    else:
        parts = np.arange(rank*psize,(rank+1)*psize)
        for idx in range(1, len(Obj)):
            NEW = Obj[idx].get_scaled_positions()
            diff = NEW[parts] - OLD[parts]
            UNWRAP_TRJ[idx][parts] = UNWRAP_TRJ[idx-1][parts] + (np.round(-diff)+diff)@LAT 
            OLD = Obj[idx].get_scaled_positions()
    UNWRAP_TRJ = comm.reduce(UNWRAP_TRJ, op=MPI.SUM, root=0)
    UNWRAP_TRJ = comm.bcast(UNWRAP_TRJ, root=0)
    return UNWRAP_TRJ

def plot_msd(MSD, TIMESTEP = 2, unit='fs', savefile=MSD_FILE):
    # TIMESTEP: femtoseconds
    fig,ax = plt.subplots()
    if unit=='fs':
        ax.plot(np.arange(0,len(MSD))*TIMESTEP,MSD)
        ax.set_xlabel("Time (fs)",fontsize=fs)
        np.savetxt(f"{unit}_"+savefile, [np.arange(0,len(MSD))*TIMESTEP, MSD], fmt='%.4f')
    elif unit=='ps':
        ax.plot(np.arange(0,len(MSD))*TIMESTEP/1000,MSD)
        ax.set_xlabel("Time (ps)",fontsize=fs)
        np.savetxt(f"{unit}_"+savefile, [np.arange(0,len(MSD))*TIMESTEP/1000,MSD], fmt='%.4f')
    ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$/fs)",fontsize=fs)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.savefig('MSD.png')
    return 0

###
if __name__ == "__main__":
    main()
