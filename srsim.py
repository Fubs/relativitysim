import sys
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
AllMovingObjects = []
AllStationaryObjects = []
class Object:
    def __init__(self, position, velocity, Color, Name):
        self.pos = np.array(position)
        self.vel = np.array(velocity)
        self.g = 1/(math.sqrt(1- self.vel[0]**2 - self.vel[1]**2))
        self.color = Color
        self.name = Name
        if velocity[0] == 0 and velocity[1] == 0:
            AllStationaryObjects.append(self)
        else:
            AllMovingObjects.append(self)

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
        # this prevents python from complaining about dividing by zero...
        if (rvx**2 + rvy**2) == 1:
            return 99999999999999999999999999999999999999999999999999
        return 1/(math.sqrt(1-(rvx**2 + rvy**2)))


# for making boosted worldlines darker than corresponding unboosted worldline, but sameish color
def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


# Some stationary objects
origin = Object([0,0,0],[0,0],"black","origin")
loc_1 = Object([10,10,0],[0,0],"red","loc_1")
loc_2 = Object([-10,10,0],[0,0],"red","loc_2")
loc_3 = Object([-10,-10,0],[0,0],"red","loc_3")
loc_4 = Object([10,-10,0],[0,0],"red","loc_4")
loc_5 = Object([0,10,0],[0,0],"red","loc_5")
loc_6 = Object([0,-10,0],[0,0],"red","loc_6")
loc_7 = Object([10,0,0],[0,0],"red","loc_7")
loc_8 = Object([-10,0,0],[0,0],"red","loc_8")

#some spaceships
ship_1 = Object([-25,0,0],[0.99,0],"cyan","ship_1")
ship_2 = Object([10,0,0],[-0.71106,0.31078],"lime","ship_2")
ship_3 = Object([20,0,0],[-0.3,0],"orange","ship_3")
ship_4 = Object([0,10,0],[0,-0.4],"violet","ship_4")



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

def generateWL(pos,vel):
    tsamples = np.linspace(0,50,50)
    xinit, yinit = pos[0], pos[1]
    vx, vy = vel[0], vel[1]
    xlist = [(xinit + vx * k) for k in tsamples]
    ylist = [(yinit + vy * k) for k in tsamples]
    return xlist, ylist, tsamples

def boostWL(X,Y,T,vel):
    stackedpts = np.stack((X,Y,T))
    boostedpts = np.dot(LorentzTransform(vel),stackedpts)
    return boostedpts[0], boostedpts[1], boostedpts[2]

def interpolateWL(iX, iY, iT, fX, fY, fT):
    numInterps = 50
    XList, YList, TList = [], [], []
    for i in range(len(iX)):
        XPath = np.linspace(iX[i],fX[i],numInterps)
        YPath = np.linspace(iY[i],fY[i],numInterps)
        TPath = np.linspace(iT[i],fT[i],numInterps)
        thisX, thisY, thisT = [], [], []
        for k in range(numInterps):
            thisX.append(XPath[k])
            thisY.append(YPath[k])
            thisT.append(TPath[k])
        XList.append(thisX)
        YList.append(thisY)
        TList.append(thisT)
    return np.array(XList).T, np.array(YList).T, np.array(TList).T



#makes the animation frames for worldlines from the perspective of viewingObject
def MakeWLS(ALLOBJECTS, viewingObject, *argv):
    otherObjects = ALLOBJECTS.copy()
    del otherObjects[ALLOBJECTS.index(viewingObject)]
    originX = viewingObject.pos[0]
    originY = viewingObject.pos[1]
    LineList = []
    for obj in otherObjects:
        shiftedX = obj.pos[0] - originX
        shiftedY = obj.pos[1] - originY
        preWLX, preWLY, preWLT = generateWL([shiftedX,shiftedY], obj.vel)
        postWLX, postWLY, postWLT = boostWL(preWLX, preWLY, preWLT, viewingObject.vel)
        XPaths, YPaths, TPaths = interpolateWL(preWLX, preWLY, preWLT, postWLX, postWLY, postWLT)
        LineList.append([XPaths,YPaths,TPaths,obj.color, obj.name])
    #now make WL for viewing object itself
    preWLX, preWLY, preWLT = generateWL([0,0], viewingObject.vel)
    postWLX, postWLY, postWLT = generateWL([0,0],[0,0])
    XPaths, YPaths, TPaths = interpolateWL(preWLX, preWLY, preWLT, postWLX, postWLY, postWLT)
    LineList.append([XPaths,YPaths,TPaths,viewingObject.color, viewingObject.name])
    titleDataString = ", original position (" + str(viewingObject.pos[0]) + ", " + str(viewingObject.pos[1]) + "), velocity (" + str(viewingObject.vel[0]) + ", " + str(viewingObject.vel[1]) + ")"
    fig.suptitle("Worldlines from perspective of " + viewingObject.name + titleDataString)
    return LineList 


    
StuffYouWannaPlot = AllMovingObjects
# LineInfo is indexed by [object][x=0, y=1, t=2][list of interp points for x,y,t]

interpedFrames=50
def PlayAnimation(perspective):
    LineInfo = MakeWLS(StuffYouWannaPlot, perspective)
    while plt.waitforbuttonpress() != True:
        continue
    for f in range(interpedFrames):
        plt.cla()
        for obj in LineInfo:
            ax.plot(obj[0][f],obj[1][f],obj[2][f], color=obj[3],label=obj[4])
        ax.set_xlim(-20,20)
        ax.set_ylim(-20,20)
        ax.set_zlim(-20,20)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("T")
        ax.legend()
        if f==0:
            plt.pause(1)
        else:
            plt.pause(0.001)

def press(event):
    print(event.key)
    sys.stdout.flush()
    if event.key == '1':
        PlayAnimation(ship_1)
    if event.key == '2':
        PlayAnimation(ship_2)
    if event.key == '3':
        PlayAnimation(ship_3)
    if event.key == '4':
        PlayAnimation(ship_4)
    if event.key == 'q':
        quit()



# cid=connection id, listens for key presses
cid = plt.gcf().canvas.mpl_connect('key_press_event', press)

#event loop
PlayAnimation(ship_1)
while True:
    plt.waitforbuttonpress()
    






