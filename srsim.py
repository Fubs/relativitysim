#from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.animation import FuncAnimation
import numpy as np
import math

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

def unitize(vec):
    v = np.array(vec)
    return v / (np.linalg.norm(v))

#Object class: self.pos() returns coords as list [x,y,t], self.vel() returns velocity vector [vx,vy]
#Velocity vector is what a stationary observer at origin would measure, and is relative to c (so [0.5,0] is 0.5c in x dir.)
#OBJECTS is list of all objects
OBJECTS = []
class Object:
    def __init__(self, position, velocity):
        self.pos = position
        self.vel = velocity
        OBJECTS.append(self)
        self.g = 1/(math.sqrt(1- self.vel[0]**2 - self.vel[1]**2))

    # relative velocity with some other object
    def relv(self, otherobj):
        vs = np.array(self.vel)
        vo = np.array(otherobj.vel)
        term1 = 1/(self.g * (1 - np.dot(vs,vo)))
        if np.dot(vs,vs) == 0:
            term2 = (self.g - 1)
        else:
            term2 = (self.g - 1) * (( np.dot(vs,vo) /np.dot(vs,vs))-1)
        term3 = vo - vs + vs * term2
        return term1 * term3

    # relative gamma factor with some other object
    def relgam(self, otherobj):
        rvx = self.relv(otherobj)[0]
        rvy = self.relv(otherobj)[1]
        return 1/(math.sqrt(1-(rvx**2 + rvy**2)))


# Some stationary objects
origin = Object([0,0,0],[0,0])
loc_1 = Object([5,5,0],[0,0])
loc_2 = Object([-5,5,0],[0,0])
loc_3 = Object([-5,-5,0],[0,0])
loc_4 = Object([5,-5,0],[0,0])

#some spaceships
shp_1 = Object([-5,0,0],[0.7,0])
shp_2 = Object([5,0,0],[-0.7156,0.3578])
shp_3 = Object([5,0,0],[-0.7,0])
shp_4 = Object([0,5,0],[0,-0.4])



# Lorentz transform. Returns the location of target in the rest frame of rest_frame_obj
def LorentzTransform(rest_frame_obj, target):
    rv = np.array(rest_frame_obj.relv(target)) #relative velocity vector
    rvm = np.dot(rv,rv) #relative velocity magnitude squared
    gam = rest_frame_obj.relgam(target) #relative gamma factor
    print(rv)
    #constructing boost matrix
    R1 = [(1+(gam - 1)*(rv[0]**2 / rvm)),((gam-1)*(rv[0]*rv[1])/rvm), ((-1)*gam*rv[0])]
    R2 = [((gam - 1)*(rv[0] * rv[1] / rvm)),(1+(gam-1)*(rv[1]**2 / rvm)), ((-1)*gam*rv[1])]
    R3 = [((-1)*gam*rv[0]), ((-1)*gam*rv[1]), (gam)]
    LTMatrix = np.array([R1,R2,R3])
    print(LTMatrix)
    return np.dot(LTMatrix, np.array(target.pos))

print(LorentzTransform(shp_4,shp_3))
print(LorentzTransform(shp_3,shp_4))




