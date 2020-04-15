#Calculate the effect of changing a spring constant on the covariance
import numpy as np
import sys

N=454

fapo=sys.argv[1]
fadj=sys.argv[2]

fhis=fapo+"_v02/IGPS_"+fapo+"_v02_1us_com_"+fadj+".dat"
fedg0=fapo+"_v02/paths/dC_"+fapo+"_v02_1us_com_"+fadj+"_adj.dat"
fout0=fapo+"_v02/paths/dCnode_"+fapo+"_v02_1us_com_"+fadj+"_adj.dat"

h=np.loadtxt(fhis)

dC=np.zeros((N,N))
delC=np.zeros((N,N))
fedg=open(fedg0,"r")
l=fedg.readline()
while l:
    l2=l.split()
    i=int(l2[0])-1
    j=int(l2[1])-1
    dC[i,j]=dC[j,i]=float(l2[2])
    delC[i,j]=delC[j,i]=float(l2[3])
    l=fedg.readline()

dCtot=-np.sum(dC*(h.clip(max=0.)),axis=0)
delCtot=np.sum(delC,axis=0)

fout=open(fout0,"w")
for i in range(N):
    fout.write("{:3d} {:19.12e} {:19.12e}\n".format(i+1,dCtot[i],delCtot[i]))
fout.close()
