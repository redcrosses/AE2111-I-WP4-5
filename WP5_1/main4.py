import numpy as np
import matplotlib.pyplot as plt
def main4(I_xx, trapezoid, stringers_pos, chord_and_span, loads, spanwise_position, design: object):
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
    def K_s(a, b):
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
    # Column buckling critical stress
    def column_buckling(MOI, A, L):
        K = 4
        E = 72.4*10**9
        critical_stress = (K*np.pi**2*E*MOI)/(L**2*A)
        return critical_stress

    def tau_cr(ks, E, nu, t, b):
        tau_cr = (np.pi**2 * ks * E) / (12 * (1 - nu**2)) * (t / b)**2
        return tau_cr
    
    def sigma_cr(kc, E, nu, t, b):
        sigma_cr = (np.pi**2 * kc * E) / (12 * (1 - nu**2)) * (t / b)**2
        return sigma_cr

    # print(trapezoid)
    # print(stringers_pos)
    for i in range(len(chord_and_span[:,0])):
        chord = chord_and_span[i,0]
        I = I_xx[i]
        M = Mx(chord_and_span[i,1])
        y_max = abs(trapezoid[1,1]*chord_and_span[i,0])
        norm_stress =  (M * y_max)/(I)
        
        avg_stress = 0 #shear force divided by (front spar height times thickness, plus rear spar height times thickness)
        shear_stress = 1.5 * avg_stress #assumed that only the spar webs carry any shear flow due to the shear force
        
        current_trapezoid = chord * trapezoid
        current_stringers = chord * stringers_pos

        plt.plot(current_trapezoid[:,0], current_trapezoid[:,1], 'o')
        plt.plot(current_stringers[:,0], current_stringers[:,1], 'o')
        plt.ylim(-3,3)
        plt.gca().set_aspect("equal", adjustable='box')
        # plt.show()



if __name__ == "__main__":
    pass

#Curve fit: Buckling coefficient for rectangular isotropic plates under shear

