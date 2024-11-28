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
        self.trapezoid = points_intersection.run([self.frontsparlength, self.rearsparlength]) #code to fit the front and rear spars into the airfoil shape. Produces the with trapezoid points
        self.chord: float = chord
        
        #airfoil
        self.x1 = np.array([0])
        self.y1 = np.array([0])
        data = np.loadtxt("WP4.2/fx60126.dat") *self.chord
        self.x1 = data[:,0]
        self.y1 = data[:,1]
        
        self.trapezoid = self.trapezoid * self.chord
        self.width: float = self.trapezoid[2,0] - self.trapezoid[1,0] #width between the front and rear spar
        self.thickness: float = skin_thickness
        
        
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

def MOI_x(box, stringer_area: float, stringer_positions, thickness: float) -> float:
    wingbox = box.trapezoid
    centroid: float = tuple(centroid_of_quadrilateral(wingbox))
    beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/box.width)
    theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/box.width)
    a = box.width /np.cos(beta)
    b = box.width /np.cos(theta)
    
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

def MOI_y(box:object, stringer_area: float, stringer_positions, thickness: float) -> float:
    wingbox = box.trapezoid
    centroid: float = tuple(centroid_of_quadrilateral(wingbox))
    beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/box.width) #slant angle of top side
    theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/box.width) #slant angle of bottom side
    a = box.width /np.cos(beta)
    b = box.width /np.cos(theta)
    
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

def J(wingbox: list[tuple], stringer_area: float, stringer_positions: list[tuple], y: float, thickness: float, chord) -> float:
    moi_x = MOI_x(wingbox, stringer_area, stringer_positions, y, thickness, chord)
    moi_y = MOI_y(wingbox, stringer_area, stringer_positions, y, thickness, chord)
    return moi_x + moi_y

#compute deflection profiles of the wing (app.D.2)

def Mx(y): 
    #WP4.1
    return y


def T(y): 
    #WP4.1
    return y


class design():
    def __init__(self, frontsparlength, rearsparlength, stringer_area, skin_thickness):
        self.span_positions = np.linspace(0, 27.47721/2 ,100)
        self.rootchord = 5.24140
        self.tipchord = 1.57714
        self.chords_along_span = np.column_stack((np.interp(self.span_positions, [0, 27.47721/2], [self.rootchord, self.tipchord]), self.span_positions))
        print(self.chords_along_span)
        self.bending_displacement: list = []
        self.torsion: list = []
        self.moi_x_list: list = []
        self.moi_y_list: list = []
        self.boxes = []

        for chord_at_span in self.chords_along_span:
            box: object = WingBox(frontsparlength,rearsparlength, chord_at_span[0], skin_thickness)
            box.makestringers(30,0.95)
            moi_x: float = MOI_x(box, stringer_area, box.stringers, box.thickness)
            moi_y: float = MOI_y(box, stringer_area, box.stringers, box.thickness)
            j = moi_x + moi_y
            def dv_dy(y: float): 
                integral, _ = quad(lambda x: -Mx(x) / E * moi_x, 0, y)
                #Mx() needs defining
                return integral
            def dtheta_dy(y: float):
                #T() needs defining
                return T(y) / (G * j)
            def v(y):
                integral, _ = quad(lambda x: dv_dy(x), 0, y)
                return integral
            def theta(y):
                integral, _ = quad(lambda x: dtheta_dy(x), 0, y)
                return integral
            self.bending_displacement.append(v(chord_at_span[1]))
            self.torsion.append(theta(chord_at_span[1]))
            self.moi_x_list.append(moi_x)
            self.moi_y_list.append(moi_y)
            if chord_at_span[0] == self.chords_along_span[0,0] or chord_at_span[0] == self.chords_along_span[-1,0]:
                self.boxes.append(box)
        # print(self.boxes[])
            
    def graph(self):
            fig1 = plt.figure()
            # print(self.boxes[0].trapezoid)
            # print(self.boxes[1].trapezoid)
            trapezoids_sized = np.vstack([self.boxes[0].trapezoid, self.boxes[1].trapezoid])
            # print(trapezoids_sized)
            #vvv first subplot vvv
            ax = fig1.add_subplot(1,2,1, projection='3d')
            ax.set(xlim3d=[0,15], ylim3d=[0,3], zlim3d=[0,30], box_aspect=(3, 3/5, 6))
            self.sweep = np.radians(25)
            # #airfoil
            airfoil = np.empty([1,2])
            data = np.loadtxt("WP4.2/fx60126.dat")
            airfoil_x = data[:,0]
            airfoil_y = data[:,1]
            airfoil_z = np.zeros(airfoil_x.shape)

            airfoil_x = np.vstack([airfoil_x*self.rootchord, airfoil_x*self.tipchord + 27.47721*np.sin(self.sweep)])
            airfoil_y = np.vstack([airfoil_y*self.rootchord, airfoil_y*self.tipchord])
            airfoil_z = np.vstack([airfoil_z, np.full_like(airfoil_z, 27.47721)])

            x,y = trapezoids_sized[:4,0], trapezoids_sized[:4,1]
            z = np.zeros(x.shape) #every four points in trapezoid_sized is one trapezoid box

            x = np.vstack([x, trapezoids_sized[4:,0] + np.sin(self.sweep)*27.47721]) 
            y = np.vstack([y, trapezoids_sized[4:,1]]) #close enough for now but it's wrong lol
            z = np.vstack([z, np.full_like(z,27.47721)])
            print(x,y,x.shape)

            airfoil = ax.plot_surface(airfoil_x, airfoil_y, airfoil_z, linewidth=0, alpha=0.25)
            surface = ax.plot_surface(x, y, z, alpha=1, linewidth=0,antialiased=False)

            # #vvv second subplot vvv idk deflections or other graphs
            # fig2 = plt.figure()
            # axs = fig2.add_subplot(2, 2, figsize=(15, 12))
            # axs = fig2.add_subplot(1,2,2)
            
            # # MOI
            # axs[0, 0].plot(self.span_positions, self.moi_x_list, label="Second moment of inertia", color='blue')
            # axs[0, 1].plot(self.span_positions, self.moi_y_list, label="Second moment of inertia", color='red')
            # axs[0, 0].set_title("Spanwise second moment of inertia in x-axis")
            # axs[0, 1].set_title("Spanwise second moment of inertia in x-axis")
            # axs[0, 0].set_ylabel("MOI_xx (m^4)")
            # axs[0, 1].set_ylabel("MOI_yy (m^4)")

            # # Displacements
            # axs[1, 0].plot(self.span_positions, self.bending_displacement, label="Bending displacement", color='blue')
            # axs[1, 1].plot(self.span_positions, self.torsion, label="Torsional twist", color='red')
            # axs[1, 0].set_title("Spanwise variation of bending displacement")
            # axs[1, 1].set_title("Spanwise variation of torsional twist")
            # axs[1, 0].set_ylabel("Bending displacement (m)")
            # axs[1, 1].set_ylabel("Torsional twist (rad)")

            # for ax in axs.flat:
            #     ax.set_xlabel("Spanwise Position (m)")
            #     ax.legend()
            #     ax.grid()

            fig1.tight_layout()
            plt.show()

# print(box.trapezoid)
design = design(0.11, 0.09, 0, 0.001)
design.graph()

#box.trapezoid provides the trapezoid points, 
#box.frontsparlength provides the front spar length
#box.rearsparlength provides the rear spar length
#box.width provides the length between the spars

#calculate the centroid, second moment of area
