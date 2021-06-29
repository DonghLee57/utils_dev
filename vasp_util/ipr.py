##### Inverse Participation Ratio calculation OCT.2009  #####
import sys,string
p=open(str(sys.argv[1]),'r') #PROCAR
g=open(str(sys.argv[2]),'w')
pp=p.readlines() #PROCAR total read

line=pp[1].split()
kpoint=int(line[3]); bands=int(line[7]) ; ions=int(line[-1])

data=[]
for i in range(len(pp)):
  if 'band ' in pp[i]:
    sum=0
    s_deno = 0
    s_numer= 0
    line=pp[i].split()
    energy=float(line[4])
    line=pp[i+3+ions].split()
    norm=float(line[-1])
    for j in range(ions):
      line=pp[i+j+3].split()
      sq_value1=float(line[-1])
      sq_value = sq_value1**2
      sum=sum+sq_value
    total=sum/norm**2
    data.append([energy, total, sum])

data.sort(key=lambda x:x[0])
for i in range(len(data)):
 g.write('%3.8f  %3.8f  %3.8f \n' %(data[i][0],data[i][1],data[i][2]))
