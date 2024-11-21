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
        self.trapezoid = points_intersection.run([self.frontsparlength, self.rearsparlength]) #code to fit the front and rear spars into the airfoil shape. Produces the with trapezoid points
        self.width: float = self.trapezoid[2,0] - self.trapezoid[1,0] #width between the front and rear spar

    def draw(self):
        #wingbox
        self.wingboxtodisp = self.trapezoid
        self.wingboxtodisp = np.append(self.trapezoid, [self.trapezoid[0]], axis=0)
        self.x2,self.y2 = zip(*self.wingboxtodisp)
        #airfoil
        self.x1 = np.array([0])
        self.y1 = np.array([0])
        data = np.loadtxt("WP4.2/fx60126.dat")
        self.x1 = data[:,0]
        self.y1 = data[:,1]
        
        plt.plot(self.x1,self.y1)
        plt.plot(self.x2,self.y2)
        plt.gca().set_aspect('equal')
        plt.show()
    def centroid(self):
        pass
    def secondmomentarea(self):
        pass
#draw/have the geomerty of a wing box 
box = WingBox(0.11,0.09) #frontspar LENGTH, rearspar LENGTH (the positions of the spars are calculated in code to fit into the airfoil)
box.draw()
#box.trapezoid provides the trapezoid points, 
#box.frontsparlength provides the front spar length
#box.rearsparlength provides the rear spar length
#box.length provides the length between the spars

#calculate the centroid, second moment of area

def MOI_x(wingbox: list[tuple], stringer_area: float, stringer_positions: list[tuple], y: float, thickness: float, chord) -> float:
    centroid: float = tuple(centroid_of_quadrilateral(wingbox))
    beta: float = np.arctan(abs(wingbox[3][1]-wingbox[0][1])/box.width)
    theta: float = np.arctan(abs(wingbox[2][1]-wingbox[1][1])/box.width)
    a = box.width/np.cos(beta)
    b = box.width/np.cos(theta)
    
    #top side
    I_xx_1 = thickness * (a**3)*(np.sin(beta)**2)*(1/12) + (thickness*a) * (abs((a/2)*np.sin(beta))+abs(wingbox[0][1]-centroid[1]))**2
    #front spar
    I_xx_2 = (box.frontsparlength * thickness**3)*(1/12) + (thickness*box.frontsparlength) * (((wingbox[0][1]+wingbox[1][1])/2)-centroid[1])**2
    #bottom side
    I_xx_3 = thickness * (b**3)*(np.sin(theta)**2)*(1/12) + (thickness*b) * (abs((b/2)*np.sin(theta))+abs(wingbox[2][1]-centroid[1]))**2
    #rear spar
    I_xx_4 = (box.rearsparlength * thickness**3)*(1/12) + (thickness*box.rearsparlength) * (((wingbox[2][1]+wingbox[3][1])/2)-centroid[1])**2
    #stringers
    I_stringers = sum([stringer_area*(abs(pos[1])-centroid[1])**2 for pos in stringer_positions])
    return I_xx_1 + I_xx_2 + I_xx_3 + I_xx_4 + I_stringers

def MOI_y(wingbox: list[tuple], stringer_area: float, stringer_positions: list[tuple], y: float, thickness: float, chord) -> float:
    centroid: float = tuple(centroid_of_quadrilateral(wingbox))
    beta: float = np.arctan(abs(wingbox[3][1]-wingbox[0][1])/box.width)
    theta: float = np.arctan(abs(wingbox[2][1]-wingbox[1][1])/box.width)
    a = box.width/np.cos(beta)
    b = box.width/np.cos(theta)
    
    #top side
    I_yy_1 = thickness * (a**3)*(np.cos(beta)**2)*(1/12) + (thickness*a) * (wingbox[0][0]+(a/2)*np.cos(beta)-centroid[0])**2
    #front spar
    I_yy_2 = (box.frontsparlength * thickness**3)*(1/12) + (thickness*box.frontsparlength) * (wingbox[1][0]-centroid[0])**2
    #bottom side
    I_yy_3 = thickness * (b**3)*(np.cos(theta)**2)*(1/12) + (thickness*b) * (wingbox[2][0]+(b/2)*np.cos(theta)-centroid[0])**2
    #rear spar
    I_yy_4 = (box.rearsparlength * thickness**3)*(1/12) + (thickness*box.rearsparlength) * (wingbox[2][0]-centroid[0])**2
    #stringers
    I_stringers = sum([stringer_area*(abs(pos[0])-centroid[0])**2 for pos in stringer_positions])
    return I_yy_1 + I_yy_2 + I_yy_3 + I_yy_4 + I_stringers

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

