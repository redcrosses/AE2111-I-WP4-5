import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad
import points_intersection
#Constraints
#- The wing tip displacement should not exceed 15% of the total span of the wing.
#- The wing tip rotation should not exceed +/- 10°.

class WingBox():
    def __init__(self, frontsparpos, rearsparpos):
        self.frontsparpos: float = frontsparpos
        self.rearsparpos: float = rearsparpos
        self.length: float = self.rearsparpos - self.frontsparpos
        
    def draw(self):
        #airfoil
        self.x1 = np.array([0])
        self.y1 = np.array([0])
        data = np.loadtxt("WP4.2/fx60126.dat")
        self.x1 = data[:,0]
        self.y1 = data[:,1]
        
        #wingbox
        self.sparpos = points_intersection.run([self.frontsparpos, self.rearsparpos])
        # print(self.sparpos)
        self.x2,self.y2 = zip(*self.sparpos)
        
        plt.plot(self.x1,self.y1)
        plt.plot(self.x2,self.y2)
        plt.gca().set_aspect('equal')
        plt.show()
    def centroid(self):
        pass
    def secondmomentarea(self):
        pass
#draw/have the geomerty of a wing box
#calculate the centroid, second moment of area

box = WingBox(0.1,0.5) #box from 10%c to 50%c
box.draw()

def bending_displacement(MOI: float, loads: list) -> list[float]:
    return

def torsion_rotation(G: float, loads: list) -> list[float]:
    return


#compute deflection profiles of the wing (app.D.2)

def Mx(y):
    return y

def EI_xx(y):
    return y

def T(y):
    return y

def GJ(y):
    return y

def dv_dy(y):
    integral, _ = quad(lambda x: -Mx(x) / EI_xx(x), 0, y)
    return integral

def v(y):
    integral, _ = quad(lambda x: dv_dy(x), 0, y)
    return integral

def dtheta_dy(y):
    return T(y) / GJ(y)

def theta(y):
    integral, _ = quad(lambda x: dtheta_dy(x), 0, y)
    return integral

