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


if __name__ == "__main__":
    pass

#Curve fit: Buckling coefficient for rectangular isotropic plates under shear