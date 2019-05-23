#python script for OVITO
#sc: tags are considered once.
#dc: tags are considered upto the its first neighbors.
#tc: #dc: tags are considered upto the its second neighbors.

from ovito.data import *
import numpy as np

def getDIS(lat,pos1,pos2):#Position arguments should be Cartensian.
 real_dist = 9999
 for a in range(-1,2):
  for b in range(-1,2):
   for c in range(-1,2):
    trans_pos = pos2 + a*lat[0] + b*lat[1] + c*lat[2]
    dist = sum((pos1-trans_pos)**2)**0.5
    if dist < real_dist:
     real_dist = dist
     Disp = trans_pos-pos1
 return real_dist, Disp
 
def modify(frame, data):
    cols = ['c_tot[3]','c_tot[5]']
    Q8 = data.particles[cols[0]]
    Q12 = data.particles[cols[1]]
    
    cr_val = 0.20
    n_nei = 3
    cutoff = 3.2
    
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
    
    finder = CutoffNeighborFinder(cutoff, data)
    for idx in range(data.number_of_particles):
        NV_list = []
        c = 0
        for neigh in finder.find(idx):
            if tag1[neigh.index] and tag2[neigh.index]: c += 1
            else: break
        if c >= n_nei:
            colors[idx] = np.array(c_code)
    
    data.particles_.create_property('Color', data=colors)
