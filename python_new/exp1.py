import numpy as np
import matplotlib.pyplot as plt
from Jabc import Jabc, Jabc2

def s2comp(theta, phi):
    return np.sin(theta)/(1+np.cos(theta)) * np.exp(phi * 1j)


pi = np.pi
GroupA=[
    (pi/6, pi/2),
    (pi/3, pi/3),
    (pi/3, 2*pi/3),
    (pi/2, 5*pi/18),
    (pi/2, pi/2),
    (pi/2, 13*pi/18),
    (2*pi/3, pi/3),
    (2*pi/3, 2*pi/3)]

GroupB=[
    (pi/6, 7*pi/6),
    (pi/3, pi),
    (pi/3, 4*pi/3),
    (pi/2, 17*pi/18),
    (pi/2, 7*pi/6),
    (pi/2, 25*pi/18),
    (2*pi/3, pi),
    (2*pi/3, 4*pi/3)]

GroupC=[
    (pi/6, -pi/6),
    (pi/3, 5*pi/3),
    (pi/3, 0),
    (pi/2, 29*pi/18),
    (pi/2, -pi/6),
    (pi/2, pi/18),
    (2*pi/3, 5*pi/3),
    (2*pi/3, 0)]

GroupD=[(5*pi/6, 0),
        (5*pi/6, pi/2),
        (5*pi/6, pi),
        (5*pi/6, pi*3/2)]

maxX, maxY, minX, minY = 0, 0, 0, 0

for i in range(6,8):
    zs_A = [s2comp(GroupA[j][0], GroupA[j][1]) for j in range(i)]
    X = [x.real for x in zs_A]
    Y = [x.imag for x in zs_A]
    maxX = max(max(X), maxX)
    maxY = max(max(Y), maxY)
    minX = min(min(X), maxX)
    minY = min(min(Y), maxY)
    plt.scatter(X,Y, color='red')
    
    zs_B = [s2comp(GroupB[j][0], GroupB[j][1]) for j in range(i)]
    X = [x.real for x in zs_B]
    Y = [x.imag for x in zs_B]
    maxX = max(max(X), maxX)
    maxY = max(max(Y), maxY)
    minX = min(min(X), maxX)
    minY = min(min(Y), maxY)
    plt.scatter(X,Y, color='green')
    
    zs_C = [s2comp(GroupC[j][0], GroupC[j][1]) for j in range(i)]
    X = [x.real for x in zs_C]
    Y = [x.imag for x in zs_C]
    maxX = max(max(X), maxX)
    maxY = max(max(Y), maxY)
    minX = min(min(X), maxX)
    minY = min(min(Y), maxY)
    plt.scatter(X,Y, color='blue')

    zs_D0 = [s2comp(GroupD[j][0], GroupD[j][1]) for j in range(len(GroupD))]
    
    zs_An = [s2comp(GroupA[j][0], GroupA[j][1]) for j in range(i,8)]
    zs_Bn = [s2comp(GroupB[j][0], GroupB[j][1]) for j in range(i,8)]
    zs_Cn = [s2comp(GroupC[j][0], GroupC[j][1]) for j in range(i,8)]

    zs_D = zs_An + zs_Bn + zs_Cn + zs_D0
    X = [x.real for x in zs_D]
    Y = [x.imag for x in zs_D]
    maxX = max(max(X), maxX)
    maxY = max(max(Y), maxY)
    minX = min(min(X), maxX)
    minY = min(min(Y), maxY)

    plt.scatter(X,Y, color='black')
    plt.xlim([-maxX, maxX])
    plt.ylim([-maxX, maxX])
    plt.show()
    
    print("Subsystem size={}".format(i))
#    print("J(A,B,C)={}".format(Jabc(zs_A, zs_B, zs_C, zs_D)))
#    print("J(A,B,C)={} (Compressed)".format(Jabc2(zs_A, zs_B, zs_C, zs_D)))
    
