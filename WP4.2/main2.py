import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad


#Constraints
#- The wing tip displacement should not exceed 15% of the total span of the wing.
#- The wing tip rotation should not exceed +/- 10°.

class WingBox():
    def __init__(self, frontsparpos, rearsparpos, length):
        self.frontsparpos: float = frontsparpos
        self.rearsparpos: float = rearsparpos
        self.length: float = length
        
    def draw(self):
        f = open("WP4.2/fx60126.dat","r")
        print(f.read())
        self.coords = [[self.frontsparpos,0],[self.frontsparpos, 1],[self.rearsparpos, 1], [self.rearsparpos, 0.5],[self.frontsparpos, 0]]
        x,y = zip(*self.coords)
        plt.plot(x,y)
        plt.show()
    def centroid(self):
        pass
    def secondmomentarea(self):
        pass
#draw/have the geomerty of a wing box
#calculate the centroid, second moment of area

box = WingBox(2,3,1)
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

