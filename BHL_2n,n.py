import numpy as np 
import matplotlib.pyplot as plt 
from scipy.linalg import eigvals
from mpl_toolkits.mplot3d import Axes3D
import time
import csv 

#Here, we consider the {4g+2,2g+1} lattices for which a Fuchsian symmetry respecting Pfaffian orientation can be constructed

g=2
Ori=[]
for i in range(1,4*g+3):
    if i<2*g+2:
       Ori.append(-1)
    else:
       Ori.append(1)

m=4*g+2
def ad(i,m): #Defining the adjecent sites
    if i==1:
        return [m,2]
    elif i==m:
        return [m-1,1]
    else: 
        return [i-1,i+1]
    
def A(b):
   A0=np.zeros((3*m,3*m), dtype=complex)
   t,a1=2**(-1/(2*m))*np.exp(-b/2),2**(-1/m)

   for i in range(1,m+1):
      for k in range(1,4):
         for j in range(1, m+1):
               for l in range(1,4):
                  I,J=3*(i-1)+k-1,3*(j-1)+l-1 #index

                  if j==ad(i,m)[0]:
                     if k==1 and l==3: 
                           A0[I][J]=-a1*Ori[j-1]
                  elif j==ad(i,m)[1]:
                     if k==3 and l==1:
                           A0[I][J]=a1*Ori[i-1]
                  else:
                     if i==j:
                           if k==1 and l==3:
                              A0[I][J]=-a1
                              A0[J][I]=a1
                           elif k==2 and l==1:
                              A0[I][J]=-t
                              A0[J][I]=t
                           else: 
                              if k==3 and l==2:
                                 A0[I][J]=-t
                                 A0[J][I]=t

   return A0

def H(p,q,r,s,b):
   j=complex(0,1)
   d=np.exp(b)
   A0=A(b)
   H0=A0
   C=[np.exp(j*p),np.exp(j*q),np.exp(j*r),np.exp(j*s),np.exp(j*(-p+q-r+s))]
   for i in range(1,2*g+1):
      I,J=3*(i-1)+2-1,3*(i+(2*g-1))+2-1
      H0[I][J]+=-C[i-1]*d
      H0[J][I]+=np.conjugate(C[i-1])*d

   Eig=np.imag(eigvals(H0))
   E_s=[e for e in Eig if e>=0]
   E0=np.sort(E_s)
   #s=np.sum(E_s)

   return E0[0]


B=np.linspace(0,1.5,50)
O=[0 for b in B]
Mesh=25

K1=np.linspace(-np.pi,np.pi,Mesh)
K2=np.linspace(-np.pi,np.pi,Mesh)
K3=np.linspace(-np.pi,np.pi,Mesh)
K4=np.linspace(-np.pi,np.pi,Mesh)


E=[]
for b in B:
   E0=[]
   for k1 in K1:
      for k2 in K2:
         for k3 in K3:
            for k4 in K4:
               E0.append(H(k1, k2, k3, k4, b))
   
   Es=np.sort(E0)
   E.append(Es[:500])

   with open('Eigenvalues_10,5.csv', 'w', newline='') as csvfile:
     csvwriter = csv.writer(csvfile)
     for b in B:
         E0=[]
         for k1 in K1:
            for k2 in K2:
               for k3 in K3:
                  for k4 in K4:
                     E0.append(H(k1, k2, k3, k4, b))

         Es=np.sort(E0)
#         # Write the current b and its corresponding E0 to the CSV file
         csvwriter.writerow(Es[:500])
         #csvwriter.writerow([s/(Mesh)**4])

# plt.plot(B,E,'r+')
# plt.plot(B,O, 'g')
# plt.show()





