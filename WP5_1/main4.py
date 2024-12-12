


import numpy as np
def main4(I_xx, trapezoids, span_and_chord, loads, spanwise_position):
    def K_c(a,b): # curve fit for skin buckling coefficient Kc
        r=a/b
        if r >0.69 and r<=1.12:
            return -104.20542*r**4+357.54147*r**3-411.62491*r**2+165.89543*r
        elif r<=1.81:
            return -8.25872*r**4+45.66709*r**3-86.31837*r**2+60.2625*r
        elif r<=5:
            return 0.049532*r**4-0.782365*r**3+4.59504*r**2-11.99811*r+-11.99811
        return "FUCK YOU"
    def Mx(y): 
        return np.interp(y, spanwise_position, loads[0][1], 0)
    def T(y): 
        return np.interp(y, spanwise_position, loads[0][2], 0)
    def K_s(a, b):
        r = a/b
        return 136.31117 - 378.14535*r + 497.60785*r**2 - 366.68125*r**3 + 163.8237*r**4 - 45.33579*r**5 + 7.595018*r**6  - 0.7056433*r**7 + 0.02790314*r**8
    # t1 and l1 are the dimensions of the rectangle parallel to the wingbox and t2 and l2 are the dimensions for the rectangle,
    #perperndicular to the wing box, and d1 is the distance between the x axis and the centroid of the parallel rectangle, 
    #higher order terms cancel out, MADE BY VICTOR BOSS
    def stringer_sizing(length_1,length_2, thickness_1, thickness_2,d1):
        Ixx =  length_1 * thickness_1 * d1 ** 2 + (1/12) * length_2 ** 3 * thickness_2
        Area_parallel = thickness_1 * length_1
        Area_perpendicular = thickness_2 * length_2 
        Total_area_stringer = Area_parallel + Area_perpendicular
        return Ixx, Area_perpendicular, Area_parallel, Total_area_stringer

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



if __name__ == "__main__":
    pass

#Curve fit: Buckling coefficient for rectangular isotropic plates under shear

