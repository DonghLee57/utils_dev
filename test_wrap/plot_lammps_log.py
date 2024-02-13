import numpy as np
import sys
import matplotlib.pyplot as plt

with open(sys.argv[1],'r') as o: tmp = o.readlines()

tag4save = False
dat = []
nions = 0
for i in range(len(tmp)):
    line = tmp[i].split()
    if len(line) >= 6:
        if line[0] == 'Step':
            tag4save = True
            idx4ste = line.index('Step')
            idx4energy = line.index('TotEng')
            continue
        elif line[0] == 'Loop':
            tag4save = False
            continue
        if tag4save:
            dat.append([int(line[0]), float(line[idx4energy])/nions])
    if 'atoms' in line and len(line) == 2:
        nions = int(line[0])
dat = np.array(dat).T

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(dat[0], dat[1],'r')

ax.set_ylabel('Per-atom energy (eV/atom)', fontsize=12)
ax.set_xlabel('Step', fontsize=12)
plt.tight_layout()
plt.show()
