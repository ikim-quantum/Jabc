from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.axes(projection="3d")

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
        


def s2c(theta, phi):
    x = np.cos(phi) * np.sin(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(theta)
    return x, y, z


xpoints=[]
ypoints=[]
zpoints=[]
color=[]

for pt in GroupA:
    x, y, z = s2c(pt[0], pt[1])
    xpoints.append(x)
    ypoints.append(y)
    zpoints.append(z)
    color.append("r")

for pt in GroupB:
    x, y, z = s2c(pt[0], pt[1])
    xpoints.append(x)
    ypoints.append(y)
    zpoints.append(z)
    color.append("g")

for pt in GroupC:
    x, y, z = s2c(pt[0], pt[1])
    xpoints.append(x)
    ypoints.append(y)
    zpoints.append(z)
    color.append("b")

for pt in GroupD:
    x, y, z = s2c(pt[0], pt[1])
    xpoints.append(x)
    ypoints.append(y)
    zpoints.append(z)
    color.append("black")

ax.scatter3D(xpoints, ypoints, zpoints, s=100, c=color)

plt.show()
