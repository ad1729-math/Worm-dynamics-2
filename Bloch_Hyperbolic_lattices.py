#Here, we study the spectrum of the hyperbolic lattice Ising model using techniques from algebraic geometry.
#We have shown that full Kasteleyn matrix can be obtained by operations of Fuchsian group elements and thus the full matrix
#reduces to only fundamental domains.

import numpy as np 
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from mpl_toolkits.mplot3d import Axes3D
import time
import csv 

b=1
g=2
j=complex(0,1)

#[(1,2,3,4),(5,6,7,8)] where i and i+4 are identified 
Ori=[1,1,1,1,1,1,1,-1] #Orientation in clockwise direction from i to i+1 (mod 8)
#Ori=[1,-1,-1,-1]
m=4*g

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
    C=[np.exp(j*p),np.exp(j*q),np.exp(j*r),np.exp(j*s)]
    for i in range(1,2*g+1):
        I,J=3*(i-1)+2-1,3*(i+(2*g-1))+2-1
        H0[I][J]+=-C[i-1]*d
        H0[J][I]+=np.conjugate(C[i-1])*d

    Eig=np.imag(eigvals(H0))
    E_s=[np.log(e) for e in Eig if e>=0]
    # E0=np.sort(E_s)
    s=np.sum(E_s)

    return s #E0[0]

#First lets fix k3,k4 and only vary k1,k2

# K1=np.linspace(0,2*np.pi,100)
# K2=np.linspace(0,2*np.pi,100)
# k3,k4=0,np.pi

# b=0.2
# Min_eig=[]
# for k1 in K1:
#     for k2 in K2:
#         Min_eig.append([k1,k2,H(k1,k2,k3,k4,b)])

# Min_eig = np.array(Min_eig)
# x, y, f_xy = Min_eig[:, 0], Min_eig[:, 1], Min_eig[:, 2]

# # Create a 3D scatter plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.scatter(x, y, f_xy, c=f_xy, cmap='viridis', marker='o')

# # Labels
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis')
# ax.set_zlabel('f(X, Y)')
# ax.set_title('3D Scatter Plot of f(x,y)')

# plt.show()                        

#All bands

st=time.time()

B=np.linspace(0,1.5,50)

Mesh=25

# K1=np.linspace(-np.pi,np.pi,Mesh)
# K2=np.linspace(-np.pi,np.pi,Mesh)
# K3=np.linspace(-np.pi,np.pi,Mesh)
# K4=np.linspace(-np.pi,np.pi,Mesh)

# with open('Pos_free_energy.csv', 'w', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
#     for b in B:
#         s=0
#         for k1 in K1:
#             for k2 in K2:
#                 for k3 in K3:
#                     for k4 in K4:
#                         s+=H(k1, k2, k3, k4, b)
        
#         # Es=np.sort(E0)
#         # Write the current b and its corresponding E0 to the CSV file
#         # csvwriter.writerow(Es[:500])
#         csvwriter.writerow([s/(Mesh)**4])

et = time.time()
print(et - st)

def read_float_list_from_csv(filename):
    float_list = []

    try:
        with open(filename, 'r') as file:
            reader = csv.reader(file)  # Default delimiter is comma
            for line in reader:
                for value in line:
                    try:
                        float_value = float(value)
                        float_list.append(float_value)
                    except ValueError:
                        print(f"Skipping invalid value: {value}")
    except IOError:
        print(f"Unable to open file: {filename}")

    return float_list

# filename='Eigenvalues.csv'
# float_list =read_float_list_from_csv(filename)

# v=500
# n=50

# D=[]
# for i in range(n):
#     D_0=[]
#     for j in range(i*v, (i+1)*v):
#         D_0.append(float_list[j])
#     D.append(D_0)

# O=[0 for b in B]
# plt.plot(B,D,'r+')
# plt.plot(B,O,'g', linewidth='1')
# plt.show()

filename='Pos_free_energy.csv'
float_list =read_float_list_from_csv(filename)

v=1
n=50

D=[]
for i in range(n):
    D_0=[]
    for j in range(i*v, (i+1)*v):
        D_0.append(float_list[j])
    D.append(D_0[0])

print(D[0]/2)

# O=[0 for b in B]
# plt.plot(B,D,'r',linewidth='3')
# plt.plot(B,O,'g', linewidth='1')
# plt.show()
