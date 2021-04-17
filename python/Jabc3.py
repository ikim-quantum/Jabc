import numpy as np
import scipy as sp
from scipy import arange, conj, prod
from scipy.sparse.linalg import LinearOperator as Lo
from scipy.sparse.linalg import svds
from scipy.sparse import csr_matrix
import math
from sympy.physics.quantum import TensorProduct as tensor
from qutip import *
import matplotlib.pyplot as plt
import itertools as itr
import time

#Listing function
def Lis(k,n):
    L=[]
    for m in range(n):
        L.append(k)
    return L

#Reshape operations
def ResAB(V,d,N,a,b,c):
    return np.reshape(V,(d**(a+b),d**(N-a-b)))

def ResBC(V,d,N,a,b,c):
    M1=np.reshape(V,(d**(a),d**(N-a)))
    M2= np.reshape(M1.T,(d**(b+c),d**(N-b-c)))
    return M2
#inverting ABCD->BC:DA
def invResBC(M,d,N,a,b,c):
    R1=np.reshape(M,(d**(N-a),d**(a))).T
    R2=np.reshape(R1,(d**N))
    return R2

#V=vector, d=local dimension, N,a,b,c,=number of sites(total, A,B,C)
def KAB(V, d, N,a,b,c):
    CAB=ResAB(V,d,N,a,b,c)
    U, S, Vh = sp.linalg.svd(CAB)
    SlogS = np.zeros((d**(a+b),d**(N-a-b)))
    for i in range(min(d**(a+b),d**(N-a-b))):
        SlogS[i, i] = sp.special.xlogy(S[i],S[i]**2)
    K2=U@SlogS@Vh
    return np.reshape(K2,(d**N))

def KBC(V, d, N,a,b,c):
    CBC=ResBC(V,d,N,a,b,c)
    U, S, Vh = sp.linalg.svd(CBC)
    SlogS = np.zeros((d**(b+c),d**(N-c-b)))
    for i in range(min(d**(b+c),d**(N-c-b))):
        SlogS[i, i] = sp.special.xlogy(S[i],S[i]**2)
    K2=U@SlogS@Vh
    return invResBC(K2,d,N,a,b,c)

#J=i<[K_AB,K_BC]\rho_{ABC}>
def JABC(V,d,N,a,b,c):
    J=(KAB(V, d, N,a,b,c).conj().T)@KBC(V, d, N,a,b,c)-(KBC(V, d, N,a,b,c).conj().T)@KAB(V, d, N,a,b,c)
    return 1j*J


#Semion coefficients: Use permutations on z to change the choice of ABC

def v(i):
    if i==-1:
        return Qobj([[1],[0]],[[2],[1]])
    else:
        return Qobj([[0],[1]],[[2],[1]])
    
def Semion4(i):
    z=[1,np.exp(1j*2*np.pi/3),np.exp(1j*4*np.pi/3),0]
    c=1
    if np.sum(i)==0:
        for a1 in range(4):
            for b1 in range (a1+1,4):
                c=c*(z[a1]-z[b1])**(i[a1]*i[b1]/2)
        return c
    else:
        return 0

def Semion6(i):
    om=np.exp(1j*2*np.pi/3)
    r=(np.sqrt(3)+1)/(np.sqrt(3)-1)
    z=[1,om,om**2,r*np.sqrt(om),r*np.sqrt(om)**3,r*np.sqrt(om)**5]
    c=1
    if np.sum(i)==0:
        for a1 in range(6):
            for b1 in range (a1+1,6):
                c=c*(z[a1]-z[b1])**(i[a1]*i[b1]/2)
        return c
    else:
        return 0

def Semion8(i):
    r=(math.sqrt(3)+1)/(math.sqrt(3)-1)
    z=[1,1j,-1,-1j,r,1j*r,-r,-1j*r]
    #z[2], z[3], z[4], z[5], z[6] = z[4], z[5], z[2], z[6], z[3]
    c=1
    if np.sum(i)==0:
        for a1 in range(8):
            for b1 in range (a1+1,8):
                c=c*(z[a1]-z[b1])**(i[a1]*i[b1]/2)
        return c
    else:
        return 0

def Semion12(i):
    phi=(math.sqrt(5)+1)/2
    th=math.atan(1/phi)
    sin=math.sin(th)
    cos=math.cos(th)
    z=[sin/(1+cos),-sin/(1+cos),1j*cos/(1+sin),-1j*cos/(1+sin),np.exp(-1j*th),np.exp(1j*th),np.exp(1j*(np.pi-th)),np.exp(-1j*(np.pi-th)),1j*cos/(1-sin),-1j*cos/(1-sin),sin/(1-cos),-sin/(1-cos)]
  #  z[3], z[4], z[5], = z[5], z[3], z[4]
    c=1
    if np.sum(i)==0:
        for a1 in range(12):
            for b1 in range (a1+1,12):
                c=c*(z[a1]-z[b1])**(i[a1]*i[b1]/2)
        return c
    else:
        return 0

Gs4=Qobj(np.zeros((16,1)),[Lis(2,4),Lis(1,4)])
for i1,i2,i3,i4 in itr.product([-1,1],[-1,1],[-1,1],[-1,1]):
    Gs42=Gs4+Semion4([i1,i2,i3,i4])*tensor(v(i1),v(i2),v(i3),v(i4))
    Gs4=Gs42
NGs4=Gs4/np.linalg.norm(Gs4)

Gs6=Qobj(np.zeros((64,1)),[Lis(2,6),Lis(1,6)])
for i1,i2,i3,i4,i5,i6 in itr.product([-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1]):
    Gs62=Gs6+Semion6([i1,i2,i3,i4,i5,i6])*tensor(v(i1),v(i2),v(i3),v(i4),v(i5),v(i6))
    Gs6=Gs62
NGs6=Gs6/np.linalg.norm(Gs6)

Gs8=Qobj(np.zeros((2**8,1)),[Lis(2,8),Lis(1,8)])
for i1,i2,i3,i4,i5,i6,i7,i8 in itr.product([-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1]):
    Gs82=Gs8+Semion8([i1,i2,i3,i4,i5,i6,i7,i8])*tensor(v(i1),v(i2),v(i3),v(i4),v(i5),v(i6),v(i7),v(i8))
    Gs8=Gs82
NGs8=Gs8/np.linalg.norm(Gs8)

#time stamp
t1=time.time()

Gs12=Qobj(np.zeros((2**12,1)),[Lis(2,12),Lis(1,12)])
for i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12 in itr.product(Lis([-1,1],repeat=12)):
    Gs122=Gs12+Semion12([i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12])*tensor(v(i1),v(i2),v(i3),v(i4),v(i5),v(i6),v(i7),v(i8),v(i9),v(i10),v(i11),v(i12))
    Gs12=Gs122
NGs12=Gs12/np.linalg.norm(Gs12)

t2=time.time()
#print(JABC(NGs12,2,12,2,2,2))
print(JABC(NGs12,2,12,2,1,3))
t3=time.time()
print(t3-t1, t2-t1,t3-t2)
