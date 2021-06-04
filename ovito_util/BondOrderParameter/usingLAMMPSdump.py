from ovito.data import *
import numpy as np

def modify(frame, data):
    fwritedat = False
    fwriteLOC = r'D:\\' # if fwritedat == True, assign the path to save data.
    looptag = ['sc','dc','dc-isolated']
    tag = 'dc-isolated'
    cols = ['c_tot[3]','c_tot[5]']
    Q8 = data.particles[cols[0]]
    Q12 = data.particles[cols[1]]

    cr_val = 0.21
    cutoff = 3.2
    n_nei  = 3
    #c_code = [0.2, 0.6, 0.8] #R, G, B
    #colors = np.ones((data.particles.count, 3))
    #data.particles_.create_property('Color', data=colors)
    trans = np.ones(data.particles.count)
    data.particles_.create_property('Transparency', data=trans)

    #
    tag1 = []
    tag2 = []
    numC = 0
    for idx in range(data.number_of_particles):
       if Q8[idx] > cr_val:
           tag1.append(True)
       else: tag1.append(False)
       if Q12[idx] > cr_val: 
           tag2.append(True)
       else: tag2.append(False)

    if tag == looptag[0]:#Consider the atom itself
        for idx in range(data.number_of_particles):
            if tag1[idx] and tag2[idx]:
                #colors[idx] = np.array(c_code)
                trans[idx] = np.array(0)
                numC += 1
    elif tag == looptag[1]:#Consider the first neigbors of the atom
        finder = CutoffNeighborFinder(cutoff, data)
        for idx in range(data.number_of_particles):
            c = 0
            for neigh in finder.find(idx):
                if tag1[neigh.index] and tag2[neigh.index]: c += 1
            if c >= n_nei:
                #colors[idx] = np.array(c_code)
                trans[idx] = np.array(0)
                numC += 1
    elif tag == looptag[2]:#After considering the first neighbors, remove isolated crystalline atoms
        checker = [False]*data.number_of_particles
        finder = CutoffNeighborFinder(cutoff, data)
        for idx in range(data.number_of_particles):
            c = 0
            for neigh in finder.find(idx):
                if tag1[neigh.index] and tag2[neigh.index]: c += 1
            if c >= n_nei:
                checker[idx] = True
                #colors[idx] = np.array(c_code)
                trans[idx] = np.array(0)
        for idx in range(data.number_of_particles):
            test =[]
            for neigh in finder.find(idx):
                test.append(checker[neigh.index])
            if True not in test:
                trans[idx] = np.array(1)
            if trans[idx] ==0: numC += 1
                
    #data.particles_.create_property('Color', data=colors)
    data.particles_.create_property('Transparency', data=trans)
    if fwritedat:
        if frame ==0:
            output = open(fwriteLOC+'num.Crystal.txt','w')
            output.write('%d\t%d\n' (%frame,numC))
        else:
            output = open(fwriteLOC+'num.Crystal.txt','a+')
            output.write('%d\t%d\n' %(frame,numC))
