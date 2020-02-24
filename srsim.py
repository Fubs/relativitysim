from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import numpy as np
import math

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def unitize(vec):
    v = np.array(vec)
    if np.linalg.norm(v) == 0:
        return vec
    return v / (np.linalg.norm(v))

#Object class: self.pos() returns coords as list [x,y,t], self.vel() returns velocity vector [vx,vy]
#Velocity vector is what a stationary observer at origin would measure, and is relative to c (so [0.5,0] is 0.5c in x dir.)
#OBJECTS is list of all objects
OBJECTS = []
class Object:
    def __init__(self, position, velocity):
        self.pos = np.array(position)
        self.vel = np.array(velocity)
        OBJECTS.append(self)
        self.g = 1/(math.sqrt(1- self.vel[0]**2 - self.vel[1]**2))

    # relative velocity with some other object
    def relv(self, otherobj):
        vs = np.array(self.vel)
        vo = np.array(otherobj.vel)
        if np.dot(vs, vo) == np.dot(vs,vs):
            return 0
        term1 = 1/(self.g * (1 - np.dot(vs,vo)))
        term2 = (self.g - 1) * (( np.dot(vs,vo) /np.dot(vs,vs))-1)
        term3 = vo - vs + vs * term2
        return np.array(term1 * term3)

    # relative gamma factor with some other object
    def relgam(self, otherobj):
        rvx = self.relv(otherobj)[0]
        rvy = self.relv(otherobj)[1]
        if (rvx**2 + rvy**2) == 1:
            return 99999999999999999999999999999999999999999999999999
        return 1/(math.sqrt(1-(rvx**2 + rvy**2)))


# Some stationary objects
origin = Object([0,0,0],[0,0])
loc_1 = Object([5,5,0],[0,0])
loc_2 = Object([-5,5,0],[0,0])
loc_3 = Object([-5,-5,0],[0,0])
loc_4 = Object([5,-5,0],[0,0])
loc_5 = Object([0,5,0],[0,0])
loc_6 = Object([0,-5,0],[0,0])
loc_7 = Object([5,0,0],[0,0])
loc_8 = Object([-5,0,0],[0,0])

#some spaceships
shp_1 = Object([-5,0,0],[0.01,0])
shp_2 = Object([5,0,0],[-0.7156,0.3578])
shp_3 = Object([5,0,0],[-0.3,0])
shp_4 = Object([0,5,0],[0,-0.4])



# Lorentz transform. Returns the boost matrix
def LorentzTransform(vel):
    vh = unitize(vel)
    gam = 1/(math.sqrt(1- vel[0]**2 - vel[1]**2))
    #constructing boost matrix
    R1 = [(1+(gam - 1)*(vh[0]**2)),((gam-1)*(vh[0] * vh[1])), ((-1)*gam*vel[0])]
    R2 = [(gam - 1)*(vh[0] * vh[1]),(1+(gam-1)*(vh[1]**2)), ((-1)*gam*vel[1])]
    R3 = [((-1)*gam*vel[0]), ((-1)*gam*vel[1]), (gam)]
    LTMatrix = np.array([R1,R2,R3])
    return LTMatrix

def plotobjects():
    for obj in OBJECTS:
        objX.append(obj.pos[0])
        objY.append(obj.pos[1])
        objZ.append(obj.pos[2])
        boostedpos = np.dot(LorentzTransform([0.8,0]),obj.pos)
        boostedX.append(boostedpos[0])
        boostedY.append(boostedpos[1])
        boostedZ.append(boostedpos[2])

boost_matrix = LorentzTransform([0.99,0])
#r = np.linspace(0,3,100)
#phi = np.linspace(0,2*np.pi,100)
#PHI, R= np.meshgrid(phi, r)
#X = R * np.cos(PHI)
#Y = R * np.sin(PHI)
#Z = np.sqrt(X**2 + Y**2)
samplingline = np.linspace(0,20,20)
X,Y = np.meshgrid(samplingline,samplingline/100)
Z = np.ones((len(samplingline),len(samplingline)))


stackedpts = np.stack((np.ravel(X),np.ravel(Y),np.ravel(Z)))
boostedpts = np.dot(boost_matrix,stackedpts).reshape(3,len(samplingline),len(samplingline))
ax.plot_surface(boostedpts[0],boostedpts[1],boostedpts[2], rstride=1, cstride=1, linewidth=0.2, antialiased=False, alpha=1)
ax.plot_surface(boostedpts[0],boostedpts[1],boostedpts[2]*2, rstride=1, cstride=1, linewidth=0.2, antialiased=False, alpha=1)
ax.plot_surface(boostedpts[0],boostedpts[1],boostedpts[2]*3, rstride=1, cstride=1, linewidth=0.2, antialiased=False, alpha=1)
ax.plot_surface(boostedpts[0],boostedpts[1],boostedpts[2]*4, rstride=1, cstride=1, linewidth=0.2, antialiased=False, alpha=1)
objX = []
objY = []
objZ = []
boostedX = []
boostedY = []
boostedZ = []

ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
ax.set_zlim(-5,5)
ax.scatter(objX,objY,objZ,color='b')
ax.scatter(boostedX,boostedY,boostedZ,color='r')


plt.show()




