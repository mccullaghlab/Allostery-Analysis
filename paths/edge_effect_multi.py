#Calculate the effect of changing a spring constant on the covariance
import numpy as np
import sys
from numpy import sqrt
from scipy.linalg import block_diag

fapo=sys.argv[1]
fadj=sys.argv[2]

fxyz=fapo+"_v02/IGPS_"+fapo+"_v02_1us_com_average_structure.dat"
fhis=fapo+"_v02/IGPS_"+fapo+"_v02_1us_com_"+fadj+".dat"
fedg0=fapo+"_v02/paths/edge_"+fapo+"_v02_1us_com_"+fadj+"_adj"
fout0=fapo+"_v02/paths/dC_"+fapo+"_v02_1us_com_"+fadj+"_adj"

N=454
N3=N*3

x=np.loadtxt(fxyz)
xhat=x[:,None]-x[None,:]
x0=sqrt(np.einsum('ijk,ijk->ij',xhat,xhat))
x0=np.divide(1.,x0,out=np.zeros_like(x0), where=x0!=0)
xhat=np.einsum('ijk,ij->ijk',xhat,x0)
xhat2=np.einsum('ijk,ijl->ijkl',xhat,xhat)

h=np.loadtxt(fhis)
### clean up negative spring constants  ###
h=(h-np.diag(np.diag(h))).clip(max=0.)
h-=np.diag(np.sum(h,axis=0))
### remove bond constants beyond 12 \AA cutoff
np.multiply(h,np.zeros_like(h),out=h,where=x0<1./12.)
### ----------------------------------  ###
h3=np.einsum('ij,ijkl->ikjl',h,xhat2)
hd=block_diag(*np.sum(h3,axis=2))
h3=np.reshape(h3,(N3,N3))-hd

e,v=np.linalg.eigh(h3)
c3=np.dot(v[:,6:]/e[6:],v[:,6:].T)

src=[49,103,129,224]
snk=[303,336,430,432]

Nran=range(N)

c3=np.reshape(c3,(N,3,N,3))
c3=np.einsum('ijkl->ikjl',c3)

dCtot=np.zeros((2,N,N))

for isrc in src:
    for isnk in snk:
        print(isrc,isnk)
        Cnm=c3[isrc,isnk]
        cnm2=np.sum(Cnm**2)

        edg=np.zeros((N,N))
        fedg=open(fedg0+".{:03d}.{:03d}.dat".format(isrc+1,isnk+1),"r")
        l=fedg.readline()
        while l:
            l2=l.split()
            i=int(l2[0])-1
            j=int(l2[1])-1
            edg[i,j]=edg[j,i]=float(l2[2])
            l=fedg.readline()

        fout=open(fout0+".{:03d}.{:03d}.dat".format(isrc+1,isnk+1),"w")
        for i in Nran:
            for j in Nran[:i]:
                UC1=np.dot(xhat[i,j],c3[i,isrc]-c3[j,isrc])
                UC2=np.dot(xhat[i,j],c3[i,isnk]-c3[j,isnk])
                dc=2.*0.6**2*(-np.sum(Cnm*np.outer(UC1,UC2)))
                
                UCi=np.dot(xhat[i,j],c3[i,i]-c3[j,i])
                UCj=np.dot(xhat[i,j],c3[i,j]-c3[j,j])
                UCU=np.dot(UCi-UCj,xhat[i,j])
                C0nm=-h[i,j]*np.outer(UC1,UC2)/(1.+h[i,j]*UCU)+Cnm
                dc0=0.6**2*(np.sum(C0nm**2)-cnm2)
        
                dCtot[0,i,j]+=dc
                dCtot[1,i,j]+=dc0
                fout.write("{:3d} {:3d} {:19.12e} {:19.12e} {:19.12e} {:19.12e}\n".format(i+1,j+1,dc,dc0,-h[i,j],edg[i,j]))
        fout.close()

edg=np.zeros((N,N))
fedg=open(fedg0+".dat","r")
l=fedg.readline()
while l:
    l2=l.split()
    i=int(l2[0])-1
    j=int(l2[1])-1
    edg[i,j]=edg[j,i]=float(l2[2])
    l=fedg.readline()

fout=open(fout0+".dat","w")
for i in Nran:
    for j in Nran[:i]:
                fout.write("{:3d} {:3d} {:19.12e} {:19.12e} {:19.12e} {:19.12e}\n".format(i+1,j+1,dCtot[0,i,j],dCtot[1,i,j],-h[i,j],edg[i,j]))
fout.close()
