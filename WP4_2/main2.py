def main2(loads: tuple, span_pos: list, n_tuple: tuple, frontsparlength: float, rearsparlength: float, horizontalsparthickness: float, verticalsparthickness: float, ribspacing: float, numberofstringers:float, stringer_width=0.03, stringer_height=0.03, thickness_1=0.001, thickness_2=0.001): #loads is a tuple, where the first element is positive loads, second element is negative loads
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from scipy.integrate import cumulative_trapezoid
    from WP4_2.centroid import centroid_of_quadrilateral
    import WP4_2.points_intersection
    from alive_progress import alive_bar
    print("\n\033[1m\033[4mCurrent Design:\033[0m")
    print("Load Distribution: {}\nFront Spar Length: {:.3f} [m]\nRear Spar Length: {:.3f} [m]\nHorizontal Spar Thickness: {:.3f} [m]\nVertical Spar Thickness: {:.3f} [m]\nRib Spacing: {:.3f} [m]\nNumber of Stringers: {:.1f}\nStringer Width: {:.3f} [m]\nStringer Height: {:.3f} [m]\nThickness 1: {:.3f} [m]\nThickness 2: {:.3f} [m]\n".format(
        n_tuple,
        frontsparlength, rearsparlength,
        horizontalsparthickness, verticalsparthickness,
        ribspacing, numberofstringers,
        stringer_width, stringer_height,
        thickness_1, thickness_2
    ))

    #Constraints
    #- The wing tip displacement should not exceed 15% of the total span of the wing.
    #- The wing tip rotation should not exceed +/- 10°.

    #CONSTANTS <3
    E = 72.4 * 10**9 #elastic modulus
    G = 27 * 10**9 #shear modulus

    class WingBox():
        def __init__(self, frontsparlength, rearsparlength, chord, hspar_thickness, vspar_thickness):
            self.frontsparlength: float = frontsparlength
            self.rearsparlength: float = rearsparlength
            self.trapezoid = WP4_2.points_intersection.run([self.frontsparlength, self.rearsparlength]) #code to fit the front and rear spars into the airfoil shape. Produces the with trapezoid points
            self.init_trapezoid = self.trapezoid
            self.chord: float = chord
            self.unitcentroid = tuple(centroid_of_quadrilateral(self.trapezoid))
            
            #airfoil
            self.x1 = np.array([0])
            self.y1 = np.array([0])
            data = np.loadtxt("WP4_2/fx60126.dat") *self.chord
            self.x1 = data[:,0]
            self.y1 = data[:,1]
            
            self.trapezoid = self.trapezoid * self.chord
            self.width: float = self.trapezoid[2,0] - self.trapezoid[1,0] #width between the front and rear spar
            self.hspar_thickness: float = hspar_thickness
            self.vspar_thickness: float = vspar_thickness
            
            
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

    def MOI_x(box, stringer_area: float, stringer_positions, hspar_thickness: float, vspar_thickness: float) -> float:
        wingbox = box.trapezoid
        centroid: float = tuple(centroid_of_quadrilateral(wingbox))
        beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/box.width)
        theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/box.width)
        a = box.width /np.cos(beta)
        b = box.width /np.cos(theta)
        
        #top side
        I_xx_1 = hspar_thickness * (a**3)*(np.sin(beta)**2)*(1/12) + (hspar_thickness*a) * (abs((a/2)*np.sin(beta))+abs(wingbox[0,1]-centroid[1]))**2
        #front spar
        I_xx_2 = (box.frontsparlength * vspar_thickness**3)*(1/12) + (vspar_thickness*box.frontsparlength) * ((((wingbox[0,1]+wingbox[1,1])/2)-centroid[1]))**2
        #bottom side
        I_xx_3 = hspar_thickness * (b**3)*(np.sin(theta)**2)*(1/12) + (hspar_thickness*b) * (abs((b/2)*np.sin(theta))+abs(wingbox[2,1]-centroid[1]))**2
        #rear spar
        I_xx_4 = (box.rearsparlength * vspar_thickness**3)*(1/12) + (vspar_thickness*box.rearsparlength) * ((((wingbox[2,1]+wingbox[3,1])/2)-centroid[1]))**2
        #stringers
        I_stringers = sum([stringer_area*((pos[1]-centroid[1]))**2 for pos in stringer_positions])
        return I_xx_1 + I_xx_2 + I_xx_3 + I_xx_4 + I_stringers

    def MOI_y(box:object, stringer_area: float, stringer_positions, hspar_thickness: float, vspar_thickness: float) -> float:
        wingbox = box.trapezoid
        centroid: float = tuple(centroid_of_quadrilateral(wingbox))
        beta: float = np.arctan(abs(wingbox[3,1]-wingbox[0,1])/box.width) #slant angle of top side
        theta: float = np.arctan(abs(wingbox[2,1]-wingbox[1,1])/box.width) #slant angle of bottom side
        a = box.width /np.cos(beta)
        b = box.width /np.cos(theta)
        
        #top side
        I_yy_1 = hspar_thickness * (a**3)*(np.cos(beta)**2)*(1/12) + (hspar_thickness*a) * ((wingbox[0,0]-centroid[0])+(a/2)*np.cos(beta))**2
        #front spar
        I_yy_2 = (box.frontsparlength * vspar_thickness**3)*(1/12) + (vspar_thickness*box.frontsparlength) * ((wingbox[1,0]-centroid[0]))**2
        #bottom side
        I_yy_3 = hspar_thickness * (b**3)*(np.cos(theta)**2)*(1/12) + (hspar_thickness*b) * ((wingbox[2,0]-centroid[0])+(b/2)*np.cos(theta))**2
        #rear spar
        I_yy_4 = (box.rearsparlength * vspar_thickness**3)*(1/12) + (vspar_thickness*box.rearsparlength) * ((wingbox[2,0]-centroid[0]))**2
        #stringers
        I_stringers = sum([stringer_area*((pos[0]-centroid[0]))**2 for pos in stringer_positions])
        return I_yy_1 + I_yy_2 + I_yy_3 + I_yy_4 + I_stringers

    class design():
        def __init__(self, frontsparlength, rearsparlength, hspar_thickness, vspar_thickness, n_stringers):
            self.span_positions = np.linspace(0, 27.47721/2, 100)
            self.rootchord = 5.24140
            self.tipchord = 1.57714
            self.chords_along_span = np.column_stack((np.interp(self.span_positions, [0, 27.47721/2], [self.rootchord, self.tipchord]), self.span_positions))
            self.n_stringers = n_stringers
            self.displacements = []
            self.frontsparlength = frontsparlength
            self.rearsparlength = rearsparlength
            self.hspar_thickness = hspar_thickness
            self.vspar_thickness = vspar_thickness
            #loadings found from diagrams
            def stringer_sizing(length_1,length_2, thickness_1, thickness_2):
                # t1 and l1 are the dimensions of the rectangle parallel to the wingbox and t2 and l2 are the dimensions for the rectangle,
                #perperndicular to the wing box, and d1 is the distance between the x axis and the centroid of the parallel rectangle, 
                #higher order terms cancel out, MADE BY VICTOR BOSS
                Area_parallel = thickness_1 * length_1
                Area_perpendicular = thickness_2 * length_2 
                Total_area_stringer = Area_parallel + Area_perpendicular
                y_centroid_stringer = (Area_perpendicular*length_2/2)/Total_area_stringer
                x_centroid_stringer = (Area_parallel*length_1/2)/Total_area_stringer
                Ixx_stringer =  Area_parallel * y_centroid_stringer ** 2 + (1/12) * length_2 ** 3 * thickness_2 + Area_perpendicular * (length_2-y_centroid_stringer) ** 2
                Iyy_stringer =  Area_perpendicular * x_centroid_stringer ** 2 + (1/12) * length_1 ** 3 * thickness_1 + Area_parallel * (length_1-x_centroid_stringer) ** 2

                return Ixx_stringer, Iyy_stringer, Total_area_stringer
            
            def ribsget(rib_spacing: float):
                pos = 0
                ribs = np.empty((1,2))
                while pos<max(self.chords_along_span[:,1]):
                    # print(pos)
                    # print(ribs)
                    chord = interp_chord(pos)
                    # print([chord, pos])
                    ribs = np.vstack((ribs, [chord,pos]))
                    pos+=rib_spacing
                    
                ribs = ribs[1::,:]
                return ribs, rib_spacing
            def interp_chord(z): #span input
                return np.interp(z, self.chords_along_span[:,1], self.chords_along_span[:,0])
            
            self.Ixx_stringer, self.Iyy_stringer, self.Total_area_stringer = stringer_sizing(stringer_width, stringer_height, thickness_1, thickness_2)
            self.ribspacing = ribspacing
            self.ribs, self.a = ribsget(self.ribspacing)  #returns chord then span placement

            with alive_bar(self.span_positions.shape[0]*2, title= "\033[96m {} \033[00m".format("WP4.2:"), bar='smooth', spinner='classic') as bar:
                for i in range(len(loads)):
                    self.boxes = []
                    bending_displacement: list = []
                    torsion: list = []
                    self.moi_x_list: list = []
                    self.moi_y_list: list = []
                    self.j_list: list = []
                    print("Load for ", n_tuple[i])
                    def Mx(y): 
                        return np.interp(y, span_pos, loads[i][1], 0)

                    def T(y): 
                        return np.interp(y, span_pos, loads[i][2], 0)
                    
                    box: object = WingBox(frontsparlength, rearsparlength, 1, hspar_thickness, vspar_thickness)
                    box.makestringers(self.n_stringers,0.95)
                    self.width = box.width
                    self.trapezoid = box.init_trapezoid
                    self.stringers = box.stringers
                
                    for chord_at_span in self.chords_along_span:
                        box: object = WingBox(frontsparlength, rearsparlength, chord_at_span[0], hspar_thickness, vspar_thickness)
                        # print(box.unitcentroid)
                        box.makestringers(self.n_stringers,0.95)
                        moi_x: float = MOI_x(box, self.Total_area_stringer, box.stringers, box.hspar_thickness, box.vspar_thickness) + self.Ixx_stringer*self.n_stringers
                        moi_y: float = MOI_y(box, self.Total_area_stringer, box.stringers, box.hspar_thickness, box.vspar_thickness) + self.Iyy_stringer*self.n_stringers
                        j = moi_x + moi_y + self.Iyy_stringer*self.n_stringers
                        self.moi_x_list.append(moi_x)
                        self.moi_y_list.append(moi_y)
                        self.j_list.append(j)
                        if chord_at_span[0] == self.chords_along_span[0,0] or chord_at_span[0] == self.chords_along_span[-1,0]:
                            self.boxes.append(box)
                        bar()

                    def fmoi_x(y):
                        return np.interp(y, self.span_positions, self.moi_x_list)
                    def fmoi_y(y):
                        return np.interp(y,self.span_positions, self.moi_y_list)
                    def fj(y):
                        return np.interp(y,self.span_positions,self.j_list)
                    
                    def dv_dy(y):
                        integrand = -Mx(y) / (E * fmoi_x(y))
                        return cumulative_trapezoid(integrand, y, initial=0) #we assume cond at wall -- no displacement slope
                    def v(y):
                        dv_dy_vals = dv_dy(y)
                        return cumulative_trapezoid(dv_dy_vals, y, initial=0) #we assume cond at wall -- no displacement
                    def dtheta_dy(y):
                        integrand =  T(y) / (G * fj(y))
                        return cumulative_trapezoid(integrand, y, initial = 0) # we assume cond at wall -- no twist slope
                    def theta(y):
                        dtheta_dy_vals = dtheta_dy(y)
                        return cumulative_trapezoid(dtheta_dy_vals, y, initial=0) # we assume cond at wall -- no twist
                    #common cum trapezoid W

                    bending_displacement = (v(self.span_positions))
                    torsion = (theta(self.span_positions))

                    self.displacements.append([bending_displacement, torsion]) #displacement first, torsion second
                    self.disp_req = 0.1*27.4277
                    self.twist_req = np.radians(10)
                    
            # print(self.boxes[])
            # print(self.displacements)
        def graph(self):
                fig1 = plt.figure(figsize=(7, 10))
                trapezoids_sized = np.vstack([np.append(self.boxes[0].trapezoid, [self.boxes[0].trapezoid[0]], axis=0), 
                                            np.append(self.boxes[1].trapezoid, [self.boxes[1].trapezoid[0]], axis=0)])
                #vvv first subplot vvv
                ax3d = fig1.add_subplot(2, 1, 1, projection='3d')  # First column
                ax3d.set(xlim3d=[0, 10], ylim3d=[0, 3], zlim3d=[0, 15], box_aspect=(2, 3/5, 3))

                self.sweep = np.radians(25)
                # #airfoil
                data = np.loadtxt("WP4_2/fx60126.dat")
                airfoil_x = data[:,0]
                airfoil_y = data[:,1]
                airfoil_z = np.zeros(airfoil_x.shape)

                ax2d = fig1.add_subplot(2,1,2)
                ax2d.set_aspect("equal")
                ax2d.plot(airfoil_x*self.rootchord, airfoil_y*self.rootchord)
                ax2d.plot(trapezoids_sized[:5,0], trapezoids_sized[:5,1])
                ax2d.plot(self.stringers[:,0]*self.rootchord, self.stringers[:,1]*self.rootchord,"ro")

                for rib in self.ribs:
                    ax3d.plot(airfoil_x*rib[0] + np.sin(self.sweep)*rib[1], airfoil_y*rib[0], np.full_like(airfoil_x, rib[1]), "r", antialiased=True)
                
                for i in range(len(self.stringers[:,0])):
                    stringer = self.stringers[i,:]
                    ax3d.plot([stringer[0]*self.rootchord, stringer[0]*self.tipchord + np.sin(self.sweep)*27.47721/2],
                              [stringer[1]*self.rootchord, stringer[1]*self.tipchord], 
                              [0, 27.47721/2], "r")

                airfoil_x = np.vstack([airfoil_x*self.rootchord, airfoil_x*self.tipchord + 27.47721/2*np.sin(self.sweep)])
                airfoil_y = np.vstack([airfoil_y*self.rootchord, airfoil_y*self.tipchord])
                airfoil_z = np.vstack([airfoil_z, np.full_like(airfoil_z, 27.47721/2)])

                x,y = trapezoids_sized[:5,0], trapezoids_sized[:5,1]
                z = np.zeros(x.shape) #every four points in trapezoid_sized is one trapezoid box

                x = np.vstack([x, trapezoids_sized[5:,0] + np.sin(self.sweep)*27.47721/2]) 
                y = np.vstack([y, trapezoids_sized[5:,1]])
                z = np.vstack([z, np.full_like(z,27.47721/2)])
                # print(x,y,x.shape)

                airfoil = ax3d.plot_surface(airfoil_x, airfoil_y, airfoil_z, linewidth=0, alpha=0.25)
                surface = ax3d.plot_surface(x, y, z, alpha=0.5, linewidth=0)

                fig1.tight_layout()
                # plt.show(block = False)
                
                # #vvv second subplot vvv
                fig2 = plt.figure(figsize=(15,10))
                gs = GridSpec(2, 3, figure=fig2, width_ratios=[1, 2, 2], height_ratios=[1, 1])            # print(self.boxes[0].trapezoid)
                # Spanwise second moment of inertia (2D plots)
                ax_moi_x = fig2.add_subplot(gs[0, 0])  
                ax_moi_x.plot(self.span_positions, self.moi_x_list, label="MOI_x", color='blue')
                ax_moi_x.set_title("Spanwise MOI_x")
                ax_moi_x.set_xlabel("Spanwise position (m)")
                ax_moi_x.set_ylabel("MOI_xx (m^4)")
                ax_moi_x.legend()
                ax_moi_x.grid()

                ax_moi_y = fig2.add_subplot(gs[1, 0])  
                ax_moi_y.plot(self.span_positions, self.j_list, label="J", color='red')
                ax_moi_y.set_title("Spanwise J")
                ax_moi_y.set_xlabel("Spanwise position (m)")
                ax_moi_y.set_ylabel("J (m^4)")
                ax_moi_y.legend()
                ax_moi_y.grid()
                for i in range(len(self.displacements)):
                    # Bending displacement pos
                    ax_bending = fig2.add_subplot(gs[0, i+1])  
                    ax_bending.plot(self.span_positions, self.displacements[i][0], label="Bending for n="+str(n_tuple[i]), color='blue')
                    ax_bending.plot(self.span_positions, np.sign(n_tuple[i])*self.disp_req*np.ones_like(self.span_positions), '--r', label='Displacement requirement')
                    ax_bending.set_title("Bending Displacement for n="+str(n_tuple[i]))
                    ax_bending.set_xlabel("Spanwise position (m)")
                    ax_bending.set_ylabel("Displacement (m)")
                    ax_bending.set_ylim(-5,5)
                    ax_bending.set_aspect("equal")
                    ax_bending.legend()
                    ax_bending.grid()

                    # Torsional twist pos
                    ax_torsion = fig2.add_subplot(gs[1, i+1])  
                    ax_torsion.plot(self.span_positions, self.displacements[i][1], label="Twist for n="+str(n_tuple[i]), color='red')
                    ax_torsion.plot(self.span_positions, -1*np.sign(n_tuple[i])*self.twist_req*np.ones_like(self.span_positions), '--r', label='Twist angle requirement')
                    ax_torsion.set_title("Torsional Twist for n="+str(n_tuple[i]))
                    ax_torsion.set_xlabel("Spanwise position (m)")
                    ax_torsion.set_ylabel("Twist (rad)")
                    ax_torsion.legend()
                    ax_torsion.grid()

                fig2.tight_layout()
                # plt.show(block=False)

        def max(self):
            print("\n\033[1m\033[4mBending/Torsion Requirements:\033[0m")
            print("\033[1mDisplacement Requirement\033[0m: +/-", self.disp_req, "[m]")
            print("Positive bending: ", self.displacements[0][0][-1], "[m]", end="")
            if self.displacements[0][0][-1] < self.disp_req:
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            print("Negative bending: ", self.displacements[1][0][-1], "[m]", end="")
            if self.displacements[1][0][-1] > -1*self.disp_req:
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            
            print("\n\033[1mTwist Requirement\033[0m: -/+", self.twist_req, "[rad]")
            print("Positive torsion: ", self.displacements[0][1][-1], "[rad]", end="")
            if self.displacements[0][1][-1] > -1*self.twist_req:
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            print("Negative torsion: ", self.displacements[1][1][-1], "[rad]", end="")
            if self.displacements[1][1][-1] < 1*self.twist_req:
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            # print("\n")
    design = design(frontsparlength, rearsparlength, horizontalsparthickness, verticalsparthickness, numberofstringers) #front spar length, rear spar length, horizontal spar thickness, vertical spar thickness, number of stringers
    design.max()
    # design.moi_x_list, design.trapezoid, design.stringers, design.chords_along_span
    return design

if __name__ == "__main__":
    
    pass

    #box.trapezoid provides the trapezoid points, 
    #box.frontsparlength provides the front spar length
    #box.rearsparlength provides the rear spar length
    #box.width provides the length between the spars

    #calculate the centroid, second moment of area
