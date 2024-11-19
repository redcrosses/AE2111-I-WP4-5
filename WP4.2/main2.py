import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad
from centroid import centroid_of_quadrilateral
import points_intersection

#Constraints
#- The wing tip displacement should not exceed 15% of the total span of the wing.
#- The wing tip rotation should not exceed +/- 10°.

class WingBox():
    def __init__(self, frontsparlength, rearsparlength):
        self.frontsparlength: float = frontsparlength
        self.rearsparlength: float = rearsparlength
        
    def draw(self):
        #airfoil
        self.x1 = np.array([0])
        self.y1 = np.array([0])
        data = np.loadtxt("WP4.2/fx60126.dat")
        self.x1 = data[:,0]
        self.y1 = data[:,1]
        
        #wingbox
        self.sparpos = points_intersection.run([self.frontsparlength, self.rearsparlength])
        self.sparpos.append(self.sparpos[0])
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

def MOI_x(wingbox: list[tuple], stringer_area: float, stringer_positions: list[tuple], y: float, thickness: float, chord) -> float:
    centroid: float = tuple(centroid_of_quadrilateral(wingbox))
    beta: float = np.arctan((wingbox[1][1]-wingbox[0][1])/box.length)
    theta: float = np.arctan((wingbox[3][1]-wingbox[4][1])/box.length)
    a = box.length/np.cos(beta)
    b = box.length/np.cos(theta)

    I_xx_1 = thickness * (a**3)*(np.sin(beta)**2) * (1/12) + (thickness*a) * ((a/2)*np.sin(beta)+(wingbox[0][1]-centroid[1]))**2
    I_xx_2 = (box.frontsparheight * thickness**3)*(1/12) + (thickness*box.frontsparheight) * (((wingbox[0][1]+wingbox[4][1])/2)-centroid[1])**2
    I_xx_3 = thickness * (b**3)*(np.sin(theta)**2) * (1/12) + (thickness*b) * ((b/2)*np.sin(theta)+(wingbox[3][1]-centroid[1]))**2
    I_xx_4 = (box.rearsparheight * thickness**3)*(1/12) + (thickness*box.rearsparheight) * (((wingbox[2][1]+wingbox[3][1])/2)-centroid[1])**2

    I_stringers = sum([stringer_area*(abs(pos[1])-centroid[1])**2 for pos in stringer_positions])
    return I_xx_1 + I_xx_2 + I_xx_3 + I_xx_4 + I_stringers

def J(y: float) -> float:
    return

def bending_displacement(MOI: float, loads: list) -> list[float]:
    return

def torsion_rotation(G: float, loads: list) -> list[float]:
    return


#compute deflection profiles of the wing (app.D.2)

def Mx(y): 
    #WP4.1
    return y

def EI_xx(y): 
    #shiyu function * E
    return y

def T(y): 
    #WP4.1
    return y

def GJ(y): 
    #shiyu function * G
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

