# use: 
#     python make_supercell.py POSCAR_original # # # > POSCAR
##### made by Joohee_Lee #####
##### ver 1.0 (20120810) #####

import sys, string

def isNum(s) :
	try :
		float(s)
		return True
	except ValueError :
		return False

try:
    poscar_original = sys.argv[1]
    expand_a = int(sys.argv[2])
    expand_b = int(sys.argv[3])
    expand_c = int(sys.argv[4])
except:
    sys.exit()

# file open
f=open(poscar_original,'r')

# title, scale, and lattice
f.readline(); line=f.readline().split(); scale=float(line[0])
line=f.readline().split(); A1=scale*float(line[0]); A2=scale*float(line[1]); A3=scale*float(line[2])
line=f.readline().split(); B1=scale*float(line[0]); B2=scale*float(line[1]); B3=scale*float(line[2])
line=f.readline().split(); C1=scale*float(line[0]); C2=scale*float(line[1]); C3=scale*float(line[2])

# name and number, or just number 
name=0
line=f.readline()
number_list=line.split()
if isNum(number_list[0])==0:
	name=line
        number_list=f.readline().split()

# total atom number calculation 
totnum=0
for i in range(len(number_list)):
	totnum=totnum+int(number_list[i])

# Selective dynamics or D/C
DorC=f.readline()
if DorC[0]=='S':
	DorC=f.readline()

# coordination list by the line order (list of lines)
coor_list=[]
for t in range(totnum):
        coor_list.append(f.readline())

# printing, title, scale, new_lattices, name_list, new number_list, Selective dynamics 
if name: print "make supercell", name.split()      #if atom names are written there 
else: print "make supercell"
print " 1.00000000000000"
print "     %(a1)16.13F %(a2)16.13F %(a3)16.13F" % {'a1':A1*expand_a,'a2':A2*expand_a,'a3':A3*expand_a} 
print "     %(b1)16.13F %(b2)16.13F %(b3)16.13F" % {'b1':B1*expand_b,'b2':B2*expand_b,'b3':B3*expand_b} 
print "     %(c1)16.13F %(c2)16.13F %(c3)16.13F" % {'c1':C1*expand_c,'c2':C2*expand_c,'c3':C3*expand_c} 
if name:                                           #if atom names are written there
	print '   ',name,    
for i in range(len(number_list)): 
        print '   ',int(number_list[i])*expand_a*expand_b*expand_c,
print ''
print "Selective dynamics"

# print if it is Direct
if DorC[0]=='D':
	print "Direct"
	for t in range(len(coor_list)):
		for i in range(expand_a):
			for j in range(expand_b):
				for k in range(expand_c):
                	                tempX=float(coor_list[t].split()[0]); tempY=float(coor_list[t].split()[1]); tempZ=float(coor_list[t].split()[2])
                        	        print " %(tx)19.16F %(ty)19.16F %(tz)19.16F   T   T   T" % {'tx':(i+tempX)/expand_a,'ty':(j+tempY)/expand_b,'tz':(k+tempZ)/expand_c}

# else if it is Cartesian
elif DorC[0]=='C':
	print "Cartesian"
	for t in range(len(coor_list)):
		for i in range(expand_a):
			for j in range(expand_b):
				for k in range(expand_c):
                	                tempX=float(coor_list[t].split()[0]); tempY=float(coor_list[t].split()[1]); tempZ=float(coor_list[t].split()[2])
                        	        print " %(tx)19.16F %(ty)19.16F %(tz)19.16F   T   T   T" % {'tx':tempX+(i*A1)+(j*B1)+(k*C1),'ty':tempY+(i*A2)+(j*B2)+(k*C2),'tz':tempZ+(i*A3)+(j*B3)+(k*C3)}
