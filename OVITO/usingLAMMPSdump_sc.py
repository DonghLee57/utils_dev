#python script for OVITO
#sc: tags are considered once.
#dc: tags are considered upto the its first neighbors.
#tc: #dc: tags are considered upto the its second neighbors.

from ovito.data import *
import numpy as np

def modify(frame, data):
    cols = ['c_tot[3]','c_tot[5]']
    Q8 = data.particles[cols[0]]
    Q12 = data.particles[cols[1]]

    cr_val = 0.20
    cutoff = 3.2

    c_code = [0.2, 0.6, 0.8] #R, G, B
    colors = np.ones((data.particles.count, 3))
    data.particles_.create_property('Color', data=colors)
    tag1 = []
    tag2 = []
    for idx in range(data.number_of_particles):
        tag1 = Q8[idx] > cr_val
        tag2 = Q12[idx] > cr_val   
        if tag1 and tag2:
            colors[idx] = np.array(c_code)

    data.particles_.create_property('Color', data=colors)
