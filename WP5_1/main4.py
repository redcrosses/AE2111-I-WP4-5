import numpy as np
import matplotlib.pyplot as plt
from alive_progress import alive_bar
class KcOutofBounds(Exception):
    pass

def main4(I_xx, trapezoid, stringers_pos, chord_and_span, loads, spanwise_position, design: object):
    ribspacing = design.ribspacing
    def K_c(a,b): # curve fit for skin buckling coefficient Kc
        r=a/b
        if r >0.68 and r<=1.12:
            return -95.34887*r**4+332.1724*r**3-390.21686*r**2+163.70145*r
        elif r<=1.81:
            return -7.39598*r**4+41.65687*r**3-80.1233*r**2+57.06951*r
        elif r<=5:
            return 0.0346329*r**4-0.557788*r**3+3.35053*r**2-9.00218*r+16.53865
        elif r>5:
            print("\033[31m K_c: r value is out of bounds (>5). a={:.3f}, b={:.3f}. K_c={:.3f} \033[0m".format(a,b,r))
            return 0.0346329*5**4-0.557788*5**3+3.35053*5**2-9.00218*5+16.53865
        
    def Vz(y,load_case=0):
        V = np.interp(y, spanwise_position, loads[0][0], 0)
        if V == 0:
            non_zero_values = np.array(loads[load_case][0])[np.array(loads[load_case][0]) != 0]
            closest_index = np.argmin(np.abs(np.array(spanwise_position) - y))
            V = non_zero_values[np.argmin(np.abs(non_zero_values - loads[load_case][0][closest_index]))]

        return V
    def Mx(y, load_case=0):
        # Interpolate the moment along the spanwise position for the given load case
        moment = np.interp(y, spanwise_position, loads[load_case][1], left=0, right=0)
        
        # If the moment is exactly 0, replace it with the closest non-zero value
        if moment == 0:
            non_zero_values = np.array(loads[load_case][1])[np.array(loads[load_case][1]) != 0]
            closest_index = np.argmin(np.abs(np.array(spanwise_position) - y))
            moment = non_zero_values[np.argmin(np.abs(non_zero_values - loads[load_case][1][closest_index]))]
        
        # Ensure the moment meets the threshold requirement
        return moment
    def T(y, load_case=0): 
        T = np.interp(y, spanwise_position, loads[0][2], 0)
        if T == 0:
            non_zero_values = np.array(loads[load_case][2])[np.array(loads[load_case][2]) != 0]
            closest_index = np.argmin(np.abs(np.array(spanwise_position) - y))
            T = non_zero_values[np.argmin(np.abs(non_zero_values - loads[load_case][2][closest_index]))]
        return T
    
    def K_s(a, b): #a is the long side, b is the short side! clamped edges
        r = a/b
        if r <= 5:
            return 136.31117 - 378.14535*r + 497.60785*r**2 - 366.68125*r**3 + 163.8237*r**4 - 45.33579*r**5 + 7.595018*r**6  - 0.7056433*r**7 + 0.02790314*r**8
        else:
            print("\033[31m K_S: r value is out of bounds (>5). a={:.3f}, b={:.3f}. K_c={:.3f} \033[0m".format(a,b,r))
            return 136.31117 - 378.14535*5 + 497.60785*5**2 - 366.68125*5**3 + 163.8237*5**4 - 45.33579*5**5 + 7.595018*5**6  - 0.7056433*5**7 + 0.02790314*5**8

    # Shear buckling critical stress
    def critical_shear_stress(ks, t, b, E=72.4 * 10 ** 9, nu=0.33):
        tau_cr = (np.pi ** 2 * ks * E) / (12 * (1 - nu ** 2)) * (t / b) ** 2
        return tau_cr

    def enclosed_area(front, rear, width):
        area = (front + rear) / 2 * width
        return area
    
    def torsion_shear_stress(front, rear, width, y):
        tau_s = T(y)/(2*enclosed_area(front, rear, width)*design.vspar_thickness)
        return tau_s
    
    def maxshear(front, rear, y):
        averageshear = Vz(y)/((front+rear)*design.vspar_thickness)
        return 1.5*averageshear

    def left_spar_shear_stress(front, rear, width, y):
        tau_left = maxshear(front, rear, y)-torsion_shear_stress(front, rear, width, y)
        return tau_left

    def right_spar_shear_stress(front, rear, width, y):
        tau_right = maxshear(front, rear, y)+torsion_shear_stress(front, rear, width, y)
        return tau_right

    # Column buckling critical stress
    def column_buckling(MOI, A, L):
        K = 4
        E = 72.4*10**9
        critical_stress = (K*np.pi**2*E*MOI)/(L**2*A)
        return critical_stress
    
    def sigma_cr(kc, t, b, E = 72.4e9, nu=0.33):
        sigma_cr = (np.pi**2 * kc * E) / (12 * (1 - nu**2)) * (t / b)**2
        return sigma_cr
    
    tip_chord = min(chord_and_span[:,0])
    root_chord = max(chord_and_span[:,0])
    # print(stringers_pos)
    def interp_chord(z): #span input
        return np.interp(z, chord_and_span[:,1], chord_and_span[:,0])
    
    def Ixx_interp(z):
        return np.interp(z, spanwise_position, I_xx)
    
    print("\n\033[1m\033[4mBuckling Analysis\33[0m")
    print("\033[01mConsidering the trailing edge panels as the most critical ones. All stresses are magnitudes.\033[0m")
    with alive_bar(len(design.ribs[:,1])-1, title= "\033[96m {} \033[00m".format("WP5.1:"), bar='smooth', spinner='classic') as bar:
        design.safetyfactors = np.empty((1,5))
        for i in range(len(design.ribs[:,1])-1):
            rib1 = design.ribs[i,:]
            chord1 = interp_chord(rib1[1])
            stringers1 = stringers_pos * chord1
            # print(stringers1)
            panel1 = stringers1[int(len(stringers1)/2)-2:int(len(stringers1)/2), :] #we assume that the max stress occurs at trailing edge of the wingbox
            # print(stringers1)
            # print(panel1)
            panel1_len = np.linalg.norm(panel1[1,:]-panel1[0,:])
            
            rib2 = design.ribs[i+1,:]
            chord2 = interp_chord(rib2[1])
            stringers2 = stringers_pos * chord2
            panel2 = stringers2[int(len(stringers2)/2)-2:int(len(stringers2)/2), :]
            panel2_len = np.linalg.norm(panel2[1,:]-panel2[0,:])
            
            design.b = (panel1_len+panel2_len)/2 #short side of panel for buckling analysis
            # print(design.a, design.b)
            
            I = Ixx_interp(rib1[1])
            M = Mx(rib1[1])
            V = Vz(rib1[1])
            y_max = abs(trapezoid[1,1]*rib1[0])
            norm_stress =  -(M * y_max)/(I)
            
            shear_stress_front = left_spar_shear_stress(chord1*design.frontsparlength, chord1*design.rearsparlength, chord1*design.width, rib1[1])
            # maxshear(chord1*design.frontsparlength, chord1*design.rearsparlength, V) #assumed that only the spar webs carry any shear flow due to the shear force
            shear_stress_rear = right_spar_shear_stress(chord1*design.frontsparlength, chord1*design.rearsparlength, chord1*design.width, rib1[1])

            Ks_front = K_s(design.a, chord1*design.frontsparlength) #a is long side, b is short side
            Ks_rear = K_s(design.a, chord1*design.rearsparlength) #a is long side, b is short side
            Kc = K_c(design.a,design.b)

            shear_buckling_stress_front = critical_shear_stress(Ks_front, design.vspar_thickness, design.frontsparlength) #spars
            shear_buckling_stress_rear = critical_shear_stress(Ks_rear, design.vspar_thickness, design.rearsparlength) #spars
            column_buckling_stress = column_buckling(design.Ixx_stringer, design.Total_area_stringer, design.a) #top panels
            skin_buckling_stress = sigma_cr(Kc, design.hspar_thickness, design.b)
            
            print("\033[36mNormal stress:\033[0m {:e}".format(norm_stress))
            print("\033[36mShear stress (front web):\033[0m {:e}".format(shear_stress_front))
            print("\033[36mShear stress (rear web):\033[0m {:e}".format(np.abs(shear_stress_rear)))
            
            print("Shear Buckling (front web) Critical Stress: {:e}".format(shear_buckling_stress_front), end="")
            if shear_stress_front < shear_buckling_stress_front:
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            print("Shear Buckling (rear web) Critical Stress: {:e}".format(shear_buckling_stress_rear), end="")
            if np.abs(shear_stress_rear) < np.abs(shear_buckling_stress_rear):
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            print("Column Buckling Critical Stress: {:e}".format(column_buckling_stress), end="")
            if norm_stress < column_buckling_stress:
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            print("Skin Buckling Critical Stress: {:e}".format(skin_buckling_stress), end="")
            if norm_stress < skin_buckling_stress:
                print("\033[32m Pass \033[0m")
            else:
                print("\033[31m Fail \033[0m")
            
            bar()
            design.safetyfactors = np.vstack((design.safetyfactors, np.array([[rib1[1], shear_buckling_stress_front/shear_stress_front, np.abs(shear_buckling_stress_rear)/shear_stress_rear, column_buckling_stress/norm_stress, skin_buckling_stress/norm_stress]])))

        #plotting
        plt.figure(figsize=(10, 6))
        plt.axhline(1, color="red", linestyle="--", label="Threshold")
        plt.plot(design.safetyfactors[:, 0][1:], design.safetyfactors[1:, 1], label=f'Shear Buckling (front spar)')
        plt.plot(design.safetyfactors[:, 0][1:], design.safetyfactors[1:, 2], label=f'Shear Buckling (rear spar)')
        plt.plot(design.safetyfactors[:, 0][1:], design.safetyfactors[1:, 3], label=f'Column Buckling')
        plt.plot(design.safetyfactors[:, 0][1:], design.safetyfactors[1:, 4], label=f'Skin Buckling')

        plt.title("Buckling Margins of Safety")
        plt.xlabel('Spanwise Position [m]')
        plt.ylabel('Safety Factor [-]')
        plt.ylim(0,15)
        plt.legend()
        # plt.yscale('log')
        plt.grid(True)

if __name__ == "__main__":
    pass

#Curve fit: Buckling coefficient for rectangular isotropic plates under shear

