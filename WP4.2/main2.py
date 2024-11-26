import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
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
    def __init__(self, frontsparlength, rearsparlength, skin_thickness):
        self.tipchord = 1.57714
        self.rootchord = 5.2414

        self.frontsparlength: float = frontsparlength
        self.rearsparlength: float = rearsparlength
        self.thickness: float = skin_thickness
        self.trapezoid = points_intersection.run([self.frontsparlength, self.rearsparlength]) #code to fit the front and rear spars into the airfoil shape. Produces unsized trapezoid points
        self.stringerarea = 0
        self.width = self.trapezoid[2,0] - self.trapezoid[1,0] #unsized width
        self.stringers = self.makestringers(self.trapezoid, 30, 0.9) #stringers for unsized cross-section

        #you can change these
        self.span_positions = np.linspace(0, 27.47721, 1000)
        self.chords = np.interp(self.span_positions, [0, 27.47721], [self.rootchord, self.tipchord])

        self.trapezoids_sized, self.MOIs = self.run()
        
    def draw(self):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        fig.tight_layout()
        ax.set(xlim3d=[0,15], ylim3d=[0,3], zlim3d=[0,30], box_aspect=(3, 3/5, 6))
        self.sweep = np.radians(25)
        # #airfoil
        self.airfoil = np.empty([1,2])
        data = np.loadtxt("WP4.2/fx60126.dat")
        airfoil_x = data[:,0]
        airfoil_y = data[:,1]
        airfoil_z = np.zeros(airfoil_x.shape)

        airfoil_x = np.vstack([airfoil_x*self.rootchord, airfoil_x*self.tipchord + 27.47721*np.sin(self.sweep)])
        airfoil_y = np.vstack([airfoil_y*self.rootchord, airfoil_y*self.tipchord])
        airfoil_z = np.vstack([airfoil_z, np.full_like(airfoil_z, 27.47721)])

        x,y = self.trapezoids_sized[0:4,0], self.trapezoids_sized[0:4,1]
        z = np.zeros(x.shape) #every four points in trapezoid_sized is one trapezoid box

        x = np.vstack([x, self.trapezoids_sized[-4:,0] + np.sin(self.sweep)*27.47721])
        y = np.vstack([y, self.trapezoids_sized[-4:,1]]) #close enough for now but it's wrong lol
        z = np.vstack([z, np.full_like(z,27.47721)])
        print(x,y,x.shape)

        airfoil = ax.plot_surface(airfoil_x, airfoil_y, airfoil_z, linewidth=0, alpha=0.25)
        surface = ax.plot_surface(x, y, z, alpha=1, linewidth=0,antialiased=False)
        manager = plt.get_current_fig_manager()

        plt.show()

    def makestringers(self, trapezoid, n, spacing_coeff):
        stringers = np.array([[],[]]) #stringer positions array
        self.topline = np.array([list(trapezoid[0]), list(trapezoid[-1])])
        self.bottomline = trapezoid[1:3,:]
        topsiden = int(n/2)
        self.stringerspacing = self.width*spacing_coeff/topsiden
        toppos = self.topline[:,0][0] + self.width*(1-spacing_coeff)
        bottomsiden = n-topsiden
        bottompos = self.bottomline[:,0][0] + self.width*(1-spacing_coeff)
        for i in range(topsiden): #make stringers on top
            ypos = np.interp(toppos, self.topline[:,0], self.topline[:,1])
            stringers = np.append(stringers, np.array([[toppos],[ypos]]), axis=1)
            toppos += self.stringerspacing
            # print("toppos: ", toppos, "ypos", ypos)
        for i in range(bottomsiden):
            ypos = np.interp(bottompos, self.bottomline[:,0], self.bottomline[:,1])
            stringers = np.append(stringers, np.array([[bottompos],[ypos]]), axis=1)
            bottompos += self.stringerspacing
        stringers = stringers.transpose()
        return stringers
        
    
    def MOI_x(self, wingbox, stringer_area: float, stringer_positions, thickness: float, width: float) -> float:
        centroid: float = tuple(centroid_of_quadrilateral(wingbox))
        beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/width)
        theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/width)
        a = self.width /np.cos(beta)
        b = self.width /np.cos(theta)
        I_xx_1 = thickness * (a**3)*(np.sin(beta)**2)*(1/12) + (thickness*a) * (abs((a/2)*np.sin(beta))+abs(wingbox[0,1]-centroid[1]))**2 #top side
        I_xx_2 = (self.frontsparlength * thickness**3)*(1/12) + (thickness*self.frontsparlength) * ((((wingbox[0,1]+wingbox[1,1])/2)-centroid[1]))**2 #front spar
        I_xx_3 = thickness * (b**3)*(np.sin(theta)**2)*(1/12) + (thickness*b) * (abs((b/2)*np.sin(theta))+abs(wingbox[2,1]-centroid[1]))**2 #bottom side
        I_xx_4 = (self.rearsparlength * thickness**3)*(1/12) + (thickness*self.rearsparlength) * ((((wingbox[2,1]+wingbox[3,1])/2)-centroid[1]))**2 #rear spar
        I_stringers = sum([stringer_area*((pos[1]-centroid[1]))**2 for pos in stringer_positions]) #stringers
        return I_xx_1 + I_xx_2 + I_xx_3 + I_xx_4 + I_stringers

    def MOI_y(self, wingbox, stringer_area: float, stringer_positions, thickness: float, width: float) -> float:
        centroid: float = tuple(centroid_of_quadrilateral(wingbox))
        beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/width) #slant angle of top side
        theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/width) #slant angle of bottom side
        a = width /np.cos(beta)
        b = width /np.cos(theta)
        I_yy_1 = thickness * (a**3)*(np.cos(beta)**2)*(1/12) + (thickness*a) * ((wingbox[0,0]-centroid[0])+(a/2)*np.cos(beta))**2 #top side
        I_yy_2 = (self.frontsparlength * thickness**3)*(1/12) + (thickness*self.frontsparlength) * ((wingbox[1,0]-centroid[0]))**2 #front spar
        I_yy_3 = thickness * (b**3)*(np.cos(theta)**2)*(1/12) + (thickness*b) * ((wingbox[2,0]-centroid[0])+(b/2)*np.cos(theta))**2 #bottom side
        I_yy_4 = (self.rearsparlength * thickness**3)*(1/12) + (thickness*self.rearsparlength) * ((wingbox[2,0]-centroid[0]))**2 #rear spar
        I_stringers = sum([stringer_area*((pos[0]-centroid[0]))**2 for pos in stringer_positions]) #stringers
        return I_yy_1 + I_yy_2 + I_yy_3 + I_yy_4 + I_stringers
    
    def J(moi_x, moi_y) -> float:
        return moi_x + moi_y    
    
    def run(self): #create the wingbox design
        trapezoids = np.empty([1,2])
        MOIs = np.empty([1,2])
        for chord in self.chords: 
            trapezoid_i = self.trapezoid * chord
            stringers_i = self.stringers * chord
            width_i = self.width * chord
            trapezoids = np.concatenate((trapezoids, trapezoid_i), axis=0) #every four points is a new trapezoid
            MOI_x = self.MOI_x(trapezoid_i, self.stringerarea, stringers_i, self.thickness, width_i)
            MOI_y = self.MOI_y(trapezoid_i, self.stringerarea, stringers_i, self.thickness, width_i)
            MOI = np.array([[MOI_x, MOI_y]])
            MOIs = np.concatenate((MOIs, MOI), axis=0)
            #BUG: after the concatenates, for some reason an extra element is added to the empty array. For now I just remove it and move on
        trapezoids = trapezoids[1:,:]
        MOIs = MOIs[1:,:]
        # print(trapezoids.shape, MOIs.shape)

        return trapezoids, MOIs
    
#compute deflection profiles of the wing (app.D.2)

def Mx(y): 
    #WP4.1
    return y

# def EI_xx_compute(wingbox: list[tuple], stringer_area: float, stringer_positions: list[tuple], y: float, thickness: float) -> float: 
#     I_xx = MOI_x(wingbox, stringer_area, stringer_positions, y, thickness)
#     EI_xx = E * I_xx
#     return EI_xx

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
box = WingBox(0.11,0.09, 0.001) #frontspar LENGTH, rearspar LENGTH (the positions of the spars are calculated in code to fit into the airfoil), airfoil skin thickness [m]
box.draw()
print(box.trapezoids_sized.shape, box.MOIs.shape)
#box.trapezoids_sized gives the trapezoids along the provided span range as a numpy array. Every four points is a new trapezoid.
#box.MOIs gives the moment of inertias. first column is MOI_x and second column is MOI_y