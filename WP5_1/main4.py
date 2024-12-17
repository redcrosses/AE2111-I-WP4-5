import numpy as np
import matplotlib.pyplot as plt
def main4(I_xx, trapezoid, stringers_pos, chord_and_span, loads, spanwise_position, design):
    def K_c(a,b): # curve fit for skin buckling coefficient Kc
        r=a/b
        if r >0.69 and r<=1.12:
            return -104.20542*r**4+357.54147*r**3-411.62491*r**2+165.89543*r
        elif r<=1.81:
            return -8.25872*r**4+45.66709*r**3-86.31837*r**2+60.2625*r
        elif r<=5:
            return 0.049532*r**4-0.782365*r**3+4.59504*r**2-11.99811*r+-11.99811
        return "FUCK YOU"
    def V(y):
        return np.interp(y, spanwise_position, loads[0][0], 0)
    def Mx(y): 
        return np.interp(y, spanwise_position, loads[0][1], 0)
    def T(y): 
        return np.interp(y, spanwise_position, loads[0][2], 0)
    def K_s(a, b): #a is the long side, b is the short side!
        r = a/b
        return 136.31117 - 378.14535*r + 497.60785*r**2 - 366.68125*r**3 + 163.8237*r**4 - 45.33579*r**5 + 7.595018*r**6  - 0.7056433*r**7 + 0.02790314*r**8

    # Shear buckling critical stress
    def critical_shear_stress(ks, t, b, E=72.4 * 10 ** 9, nu=0.33, ):
        tau_cr = (np.pi ** 2 * ks * E) / (12 * (1 - nu ** 2)) * (t / b) ** 2
        return tau_cr

    def enclosed_area(trapezoid):
        x = trapezoid[:, 0]
        y = trapezoid[:, 1]
        area = 0.5 * abs(
            np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))
        )
        return area
    
    def torsion_shear_stress(area, y, thickness):
        tau_s = T(y)/(2*area*thickness)
        return tau_s

    def enclosed_area(trapezoid):
        x = trapezoid[:, 0]
        y = trapezoid[:, 1]
        area = 0.5 * abs(
            np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))
        )

        return area
    def torsion_shear_stress(area, y, thickness):
        tau_s = T(y)/(2*area*thickness)
        return tau_s
    
    def maxshear():
        frontlength = design.chords_along_span[:,0]*design.frontsparlength
        rearlength = design.chords_along_span[:,0]*design.rearsparlength
        averageshear = V/((frontlength+rearlength)*design.vspar_thickness)
        maxshear = 1.5*averageshear
        
    # Column buckling critical stress
    def column_buckling(MOI, A, L):
        K = 4
        E = 72.4*10**9
        critical_stress = (K*np.pi**2*E*MOI)/(L**2*A)
        return critical_stress

    def tau_cr(ks, E, nu, t, b):
        tau_cr = (np.pi**2 * ks * E) / (12 * (1 - nu**2)) * (t / b)**2
        return tau_cr
    
    def sigma_cr(kc, t, b, E = 72.4*10**9, nu=0.33):
        sigma_cr = (np.pi**2 * kc * E) / (12 * (1 - nu**2)) * (t / b)**2
        return sigma_cr
    
    tip_chord = min(chord_and_span[:,0])
    root_chord = max(chord_and_span[:,0])
    # print(stringers_pos)
    def interp_chord(z): #span input
        return np.interp(z, chord_and_span[:,1], chord_and_span[:,0])
    
    def Ixx_interp(z):
        return np.interp(z, spanwise_position, I_xx)
    
    def ribsget(rib_spacing: float):
        pos = 0
        ribs = np.empty((1,2))
        while pos<max(chord_and_span[:,1]):
            # print(pos)
            # print(ribs)
            chord = interp_chord(pos)
            # print([chord, pos])
            ribs = np.vstack((ribs, [chord,pos]))
            pos+=rib_spacing
            
        ribs = ribs[1::,:]
        return ribs, rib_spacing
    
    ribs, a = ribsget(2) #2 meter spacing; returns chord then span placement
    # print(ribs_placement)
    # print(stringers_pos)
    # print(trapezoid)
    # print(stringers_pos)
    for i in range(len(ribs[:,1])-1):
        rib1 = ribs[i,:]
        chord1 = interp_chord(rib1[1])
        stringers1 = stringers_pos * chord1
        # print(stringers1)
        panel1 = stringers1[int(len(stringers1)/2)-2:int(len(stringers1)/2), :] #we assume that the max stress occurs at trailing edge of the wingbox
        # print(stringers1)
        # print(panel1)
        panel1_len = np.linalg.norm(panel1[1,:]-panel1[0,:])
        
        rib2 = ribs[i+1,:]
        chord2 = interp_chord(rib2[1])
        stringers2 = stringers_pos * chord2
        panel2 = stringers2[int(len(stringers2)/2)-2:int(len(stringers2)/2), :]
        panel2_len = np.linalg.norm(panel2[1,:]-panel2[0,:])
        
        b = (panel1_len+panel2_len)/2 #short side of panel for buckling analysis
        print(a, b)
        
        I = Ixx_interp(rib2[1])
        M = Mx(ribs[i+1,1])
        y_max = abs(trapezoid[1,1]*ribs[i+1,0])
        norm_stress =  (M * y_max)/(I)
        
        avg_stress = 0 #shear force divided by (front spar height times thickness, plus rear spar height times thickness)
        shear_stress = 1.5 * avg_stress #assumed that only the spar webs carry any shear flow due to the shear force

        Ks = K_s(a,b) #a is long side, b is short side
        Kc = K_c(a,b)
        shear_buckling_stress = critical_shear_stress(Ks, design.hspar_thickness, b)
        column_buckling_stress = column_buckling(design.Ixx_stringer, design.Total_area_stringer, a)
        skin_buckling_stress = sigma_cr(Kc, design.hspar_thickness,b)
        

        # plt.plot(current_trapezoid[:,0], current_trapezoid[:,1], 'o')
        # plt.plot(current_stringers[:,0], current_stringers[:,1], 'o')
        # plt.ylim(-3,3)
        # plt.gca().set_aspect("equal", adjustable='box')
        # plt.show()



if __name__ == "__main__":
    pass

#Curve fit: Buckling coefficient for rectangular isotropic plates under shear

