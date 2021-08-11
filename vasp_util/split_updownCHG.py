#To split up & down spin charge density 
# python *.py PARCHG
# >> PARCHG_up, PARCHG_down 
# You can visualize these file using VESTA.
#
# Format setting may be needed for other types of files.
#
import sys
import math

parchg = sys.argv[1]

opar = open(parchg,'r')
par_up = open("PARCHG_up",'w')
par_down = open("PARCHG_down",'w')


for i in range(6):
 ll = opar.readline()
 par_up.write("%s" %ll)
 par_down.write("%s" %ll)

ll = opar.readline()
par_up.write("%s" %ll)
par_down.write("%s" %ll)
Natom_list = map(int,ll.split())
Natom = sum(Natom_list)

for i in range(Natom+2):
 ll = opar.readline()
 par_up.write("%s" %ll)
 par_down.write("%s" %ll)

ll = opar.readline()
par_up.write("%s" %ll)
par_down.write("%s" %ll)

grid_list = map(int, ll.split())
Ngrid = grid_list[0]*grid_list[1]*grid_list[2]
Nline = int(math.ceil(Ngrid/10.0))

UpD = []
UmD = []

for i in range(Nline):
 UpD.append(map(float,opar.readline().split()))

opar.readline()
opar.readline()
opar.readline()

for i in range(Nline):
 UmD.append(map(float,opar.readline().split()))

for i in range(Nline):
 if (i == Nline-1):
  for j in range(len(UmD[i])):
   if (j == len(UmD[i])):
    par_up.write("\t%f\n"%((UpD[i][j] + UmD[i][j])/2))
    par_down.write("\t%f\n"%((UpD[i][j] - UmD[i][j])/2))
   else: 
    par_up.write("\t%f"%((UpD[i][j] + UmD[i][j])/2))
    par_down.write("\t%f"%((UpD[i][j] - UmD[i][j])/2))
 else:
  for j in range(10):
   if (j == 9):
    par_up.write("\t%f\n"%((UpD[i][j] + UmD[i][j])/2))
    par_down.write("\t%f\n"%((UpD[i][j] - UmD[i][j])/2))
   else: 
    par_up.write("\t%f"%((UpD[i][j] + UmD[i][j])/2))
    par_down.write("\t%f"%((UpD[i][j] - UmD[i][j])/2))

par_up.close()
par_down.close()
