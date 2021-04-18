import numpy as np
from Jabc import Jabc, Jabc2
import matplotlib.pyplot as plt

omega = np.exp(np.pi * 1j * 2 / 3)
r = (np.sqrt(3)+1) / (np.sqrt(3)-1)

# 4 qubits
zs = [1, omega, omega**2, 0]
zs_A = [zs[0]]
zs_B = [zs[1]]
zs_C = [zs[2]]
zs_D = [zs[3]]

X = [z.real for z in zs]
Y = [z.imag for z in zs]
color= ["red", "green", "blue", "black"]
plt.scatter(X,Y, c=color)
plt.show()

print("4 qubit case. J(A,B,C)={}".format(Jabc(zs_A, zs_B, zs_C, zs_D, verbose=False)))
print("4 qubit case. (compressed) J(A,B,C)={}".format(Jabc2(zs_A, zs_B, zs_C, zs_D, verbose=False)))

# 6 qubits
zs = [1, omega, omega**2, r * np.sqrt(omega), r * np.sqrt(omega) * omega, r * np.sqrt(omega) * omega * omega]
zs_A = [zs[0]]
zs_B = [zs[1]]
zs_C = [zs[2]]
zs_D = zs[3:]

X = [z.real for z in zs]
Y = [z.imag for z in zs]
color= ["red", "green", "blue", "black", "black", "black"]
plt.scatter(X,Y, c=color)
plt.show()
print("6 qubit case. J(A,B,C)={}".format(Jabc(zs_A, zs_B, zs_C, zs_D, verbose=False)))
print("6 qubit case. (compressed) J(A,B,C)={}".format(Jabc2(zs_A, zs_B, zs_C, zs_D, verbose=False)))


# Fibonacci
N=21
zs = []
phi= (1+np.sqrt(5))/2
for i in range(N):
    t1 = i/N
    t2 = i/phi
    zs.append(np.sqrt(t1) * np.exp(2*np.pi * 1j * t2))

X = [z.real for z in zs]
Y = [z.imag for z in zs]
plt.scatter(X,Y)
plt.show()
