import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad
from centroid import centroid_of_quadrilateral
import points_intersection

#Constraints
#- The wing tip displacement should not exceed 15% of the total span of the wing.
#- The wing tip rotation should not exceed +/- 10°.

#CONSTANTS <3
E = 72.4 * 10**9
G = 27 * 10**9

class WingBox():
    def __init__(self, frontsparlength, rearsparlength, chord, skin_thickness):
        self.frontsparlength: float = frontsparlength
        self.rearsparlength: float = rearsparlength
        self.thickness: float = skin_thickness
        self.trapezoid = points_intersection.run([self.frontsparlength, self.rearsparlength]) #code to fit the front and rear spars into the airfoil shape. Produces the with trapezoid points
        # self.stringers = self.makestringers(30,0.95) #number of stringers, how much of wingbox to cover in percent

        # #airfoil
        # self.x1 = np.array([0])
        # self.y1 = np.array([0])
        # data = np.loadtxt("WP4.2/fx60126.dat") *chord
        # self.x1 = data[:,0]
        # self.y1 = data[:,1]
        
    def draw(self):
        #wingbox
        self.wingboxtodisp = self.trapezoid
        self.wingboxtodisp = np.append(self.trapezoid, [self.trapezoid[0]], axis=0)
        self.x2,self.y2 = zip(*self.wingboxtodisp)
        plt.plot(self.stringers[:,0], self.stringers[:,1], "o")
        plt.plot(self.x1,self.y1)
        plt.plot(self.x2,self.y2)
        plt.gca().set_aspect('equal')
        plt.show()
    def makestringers(self, n, spacing_coeff):
        self.stringers = np.array([[],[]]) #stringer positions array
        self.topline = np.array([list(self.trapezoid[0]), list(self.trapezoid[-1])])
        self.bottomline = self.trapezoid[1:3,:]
        topsiden = int(n/2)
        self.stringerspacing = self.width*spacing_coeff/topsiden
        toppos = self.topline[:,0][0] + self.width*(1-spacing_coeff)
        bottomsiden = n-topsiden
        bottompos = self.bottomline[:,0][0] + self.width*(1-spacing_coeff)
        for i in range(topsiden): #make stringers on top
            ypos = np.interp(toppos, self.topline[:,0], self.topline[:,1])
            self.stringers = np.append(self.stringers, np.array([[toppos],[ypos]]), axis=1)
            toppos += self.stringerspacing
            # print("toppos: ", toppos, "ypos", ypos)
        for i in range(bottomsiden):
            ypos = np.interp(bottompos, self.bottomline[:,0], self.bottomline[:,1])
            self.stringers = np.append(self.stringers, np.array([[bottompos],[ypos]]), axis=1)
            bottompos += self.stringerspacing
        self.stringers = self.stringers.transpose()
        # print(self.stringers)
        pass
    
    def MOI_x(self, stringer_area: float, stringer_positions, thickness: float) -> float:
        wingbox = self.trapezoid
        centroid: float = tuple(centroid_of_quadrilateral(wingbox))
        beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/self.width)
        theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/self.width)
        a = self.width /np.cos(beta)
        b = self.width /np.cos(theta)
        
        #top side
        I_xx_1 = thickness * (a**3)*(np.sin(beta)**2)*(1/12) + (thickness*a) * (abs((a/2)*np.sin(beta))+abs(wingbox[0,1]-centroid[1]))**2
        #front spar
        I_xx_2 = (box.frontsparlength * thickness**3)*(1/12) + (thickness*box.frontsparlength) * ((((wingbox[0,1]+wingbox[1,1])/2)-centroid[1]))**2
        #bottom side
        I_xx_3 = thickness * (b**3)*(np.sin(theta)**2)*(1/12) + (thickness*b) * (abs((b/2)*np.sin(theta))+abs(wingbox[2,1]-centroid[1]))**2
        #rear spar
        I_xx_4 = (box.rearsparlength * thickness**3)*(1/12) + (thickness*box.rearsparlength) * ((((wingbox[2,1]+wingbox[3,1])/2)-centroid[1]))**2
        #stringers
        I_stringers = sum([stringer_area*((pos[1]-centroid[1]))**2 for pos in stringer_positions])
        return I_xx_1 + I_xx_2 + I_xx_3 + I_xx_4 + I_stringers

    def MOI_y(self, stringer_area: float, stringer_positions, thickness: float) -> float:
        wingbox = self.trapezoid
        centroid: float = tuple(centroid_of_quadrilateral(wingbox))
        beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/self.width) #slant angle of top side
        theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/self.width) #slant angle of bottom side
        a = self.width /np.cos(beta)
        b = self.width /np.cos(theta)
        
        #top side
        I_yy_1 = thickness * (a**3)*(np.cos(beta)**2)*(1/12) + (thickness*a) * ((wingbox[0,0]-centroid[0])+(a/2)*np.cos(beta))**2
        #front spar
        I_yy_2 = (box.frontsparlength * thickness**3)*(1/12) + (thickness*box.frontsparlength) * ((wingbox[1,0]-centroid[0]))**2
        #bottom side
        I_yy_3 = thickness * (b**3)*(np.cos(theta)**2)*(1/12) + (thickness*b) * ((wingbox[2,0]-centroid[0])+(b/2)*np.cos(theta))**2
        #rear spar
        I_yy_4 = (box.rearsparlength * thickness**3)*(1/12) + (thickness*box.rearsparlength) * ((wingbox[2,0]-centroid[0]))**2
        #stringers
        I_stringers = sum([stringer_area*((pos[0]-centroid[0]))**2 for pos in stringer_positions])
        return I_yy_1 + I_yy_2 + I_yy_3 + I_yy_4 + I_stringers
    
    def J(moi_x, moi_y) -> float:
        return moi_x + moi_y    
    
    def run(self):
        span_positions = np.linspace(0, 27.47721 ,100)
        chords = np.interp(span_positions, [0, 27.47721], [5.24140, 1.57714])
        trapezoids_x = np.array([])
        trapezoids_y = np.array([])
        widths = np.array([])
        print(self.trapezoid[:,0])
        for chord in chords: 
            trapezoids_x = np.append(trapezoids_x, self.trapezoid[:,0] * chord)
            trapezoids_y = np.append(trapezoids_y, self.trapezoid[:,1] * chord)
            # allstringers = self.stringers 
            # MOI_xs = self.MOI_x(0, box.stringers, box.thickness)
            # MOI_ys = self.MOI_y(0, box.stringers, box.thickness)

        # plt.scatter(trapezoids_x, trapezoids_y, span_positions)
        # plt.show()
#compute deflection profiles of the wing (app.D.2)

def Mx(y): 
    #WP4.1
    return y

def EI_xx_compute(wingbox: list[tuple], stringer_area: float, stringer_positions: list[tuple], y: float, thickness: float) -> float: 
    I_xx = MOI_x(wingbox, stringer_area, stringer_positions, y, thickness)
    EI_xx = E * I_xx
    return EI_xx

def T(y): 
    #WP4.1
    return y

def GJ(wingbox: list[tuple], stringer_area: float, stringer_positions: list[tuple], y: float, thickness: float) -> float:
    J = J(wingbox, stringer_area, stringer_positions, y, thickness)
    return G * J

#deflection functions

def dv_dy(y): 
    integral, _ = quad(lambda x: -Mx(x) / EI_xx_compute(x), 0, y)
    return integral

def v(y):
    integral, _ = quad(lambda x: dv_dy(x), 0, y)
    return integral

def dtheta_dy(y):
    return T(y) / GJ(y)

def theta(y):
    integral, _ = quad(lambda x: dtheta_dy(x), 0, y)
    return integral

#draw/have the geomerty of a wing box 
box = WingBox(0.11,0.09, 2, 0.001) #frontspar LENGTH, rearspar LENGTH (the positions of the spars are calculated in code to fit into the airfoil)
box.run()
#box.trapezoid provides the trapezoid points, 
#box.frontsparlength provides the front spar length
#box.rearsparlength provides the rear spar length
#box.width provides the length between the spars

#calculate the centroid, second moment of area
