#python script for OVITO
#sc: tags are considered once.
#dc: tags are considered upto the its first neighbors.

from ovito.data import *
import numpy as np

def modify(frame, data):
    looptag = ['sc','dc']
    tag = 'dc'
    cols = ['c_tot[3]','c_tot[5]']
    Q8 = data.particles[cols[0]]
    Q12 = data.particles[cols[1]]

    cr_val = 0.20
    cutoff = 3.2
    n_nei  = 3
    c_code = [0.2, 0.6, 0.8] #R, G, B

    colors = np.ones((data.particles.count, 3))
    data.particles_.create_property('Color', data=colors)
    tag1 = []
    tag2 = []

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
                colors[idx] = np.array(c_code)
    elif tag == looptag[1]:#Consider the first neigbors of the atom
        finder = CutoffNeighborFinder(cutoff, data)
        for idx in range(data.number_of_particles):
            c = 0
            for neigh in finder.find(idx):
                if tag1[neigh.index] and tag2[neigh.index]: c += 1
                else: break
            if c >= n_nei:
                colors[idx] = np.array(c_code)

   data.particles_.create_property('Color', data=colors)
