# python 3.X
# *Complete function readPROCAR
# *Complete function readOUTCAR
#
# *LDOS in slab analysis code
#  1) load DOSCAR (atom-resolved)
#  2) load OUTCAR? (atom coordination)
#  3) bin in x,y, or z direction
#  4) plot pdos of the allocated atoms
#
import numpy as np
import matplotlib.pyplot as plt

def main():
    #test for "readPOSCAR" function
    #print(info_pos)
    lat, types, atoms, *info_pos = readPOSCAR('POSCAR')
    print(fit_density(8.14,lat, types, atoms))

    #test for "readDOSCAR" function
    #ISPIN = 2
    #mask_u = np.arange(1,19,2)
    #mask_d = np.arange(2,19,2)
    #TDOS, PDOS, fermi = readDOSCAR("DOSCAR", ISPIN)
    #
    #fig, ax = plt.subplots()
    #ax.axhline(0,c='k')
    #ax.axvline(0,c='k',ls=':')
    #ax.plot(TDOS.T[0]-fermi,TDOS.T[1],'r')
    #ax.plot(TDOS.T[0]-fermi,TDOS.T[2]*(-1),'r')
    #ax.plot(PDOS[-1].T[0]-fermi,np.sum(PDOS[-1][:].T[mask_u],axis=0),'b')
    #ax.plot(PDOS[-1].T[0]-fermi,np.sum(PDOS[-1][:].T[mask_d],axis=0)*(-1),'b')
    #ax.set_xlim([-2,2])
    #ax.set_xlabel(r'E-E$_{f}$')
    #ax.set_ylabel('a.u.')

    #plt.tight_layout()
    #plt.show()

    #info = readPROCAR(winPATH + "PROCAR_example")

    return 1

#--------------------------------------------------------------------------------------
def readPOSCAR(filename):
    def fillcoords(data, NIONS, coords_array):
        for i in range(NIONS):
            items = list(map(float,data[i].split()[:3]))
            for j in range(3):
                coords_array[i][j] = items[j]
        return coords_array

    o = open(filename,'r')
    o.readline()
    scale = float(o.readline())
    lat = np.zeros((3,3))
    for i in range(3):
        items = list(map(float,o.readline().split()))
        for j in range(3):
            lat[i][j] = items[j]
    lat = scale*lat
    ilat = np.linalg.inv(lat)

    types = o.readline().split()
    atoms = list(map(int,o.readline().split()))
    NIONS = np.sum(atoms)

    Dcoords = np.empty((NIONS,3))
    Ccoords = np.empty((NIONS,3))
    coordtype = o.readline()[0]
    if coordtype == 'S':
        coordtype = o.readline()[0]
        if coordtype == 'D':
            Dcoords = fillcoords(o.readlines()[:NIONS], NIONS, Dcoords)
            Ccoords = np.matmul(Dcoords,lat)
        elif coordtype == 'C':
            Ccoords = fillcoords(o.readlines()[:NIONS], NIONS, Ccoords)
            Dcoords = np.matmul(Ccoords,ilat)
    else:
        if coordtype == 'D':
            Dcoords = fillcoords(o.readlines()[:NIONS], NIONS, Dcoords)
            Ccoords = np.matmul(Dcoords,lat)
        elif coordtype == 'C':
            Ccoords = fillcoords(o.readlines()[:NIONS], NIONS, Ccoords)
            Dcoords = np.matmul(Ccoords,ilat)
    return lat, types, atoms, Dcoords, Ccoords

def readDOSCAR(filename, ISPIN):
    # LORBIT = 11 format
    o = open(filename,'r')
    tmp = list(map(int,o.readline().split()))
    NIONS = tmp[0]
    FLAG_pdos = tmp[-2] > 0
    o.readline(); o.readline(); o.readline(); o.readline();
    EMAX, EMIN, NEDOS, fermi, _  = list(map(float,o.readline().split()))
    NEDOS = int(NEDOS)
    orbitals = ISPIN*9

    TDOS = np.empty((NEDOS,ISPIN+1))
    for i in range(NEDOS):
        TDOS[i] = list(map(float, o.readline().split()))[:ISPIN+1]

    if FLAG_pdos:
        PDOS = np.empty((NIONS,NEDOS,orbitals+1))
        for n in range(NIONS):
            o.readline()
            for i in range(NEDOS):
                PDOS[n][i] = list(map(float, o.readline().split()))
        return TDOS, PDOS, fermi
    else:
        return TDOS, 0, fermi

def fit_density(target_rho, lat, types, atoms):
   #target: The target density (unit: g/cm3)
   Na = 6.022E+23
   dic_mass = {'O':16.00,'Hf':178.49}
   V = np.linalg.det(lat)
   Tot_mass = 0
   for i in range(len(types)):
      Tot_mass += dic_mass[types[i]]*atoms[i]
   # Unit: g/cm3
   given_rho = Tot_mass/V/Na*10**24
   fit_lat = (target_rho/given_rho)**(1/3.)*lat
   return fit_lat

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
