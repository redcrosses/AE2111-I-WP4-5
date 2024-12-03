import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy as sp
from scipy import interpolate
# Constants
V = 10 # [m/s]
q = 0.5 * 1.225 * V**2
wing_span = 27.47721  # Total span of the wing (m)
alpha_a = 14 * np.pi / 180  # Angle of attack in radians
rho = 1.225  # Air density in kg/m^3 (sea level standard)
velocity = 0.85 * 343  # Freestream velocity in m/s

# Read AVL Data
def read_xflr_data(file):
    data_start = False
    data = []
    with open(file, 'r') as f:
        for line in f:
            if "Main Wing" in line and not data_start:
                data_start = True
                continue
            if data_start and (not line.strip() or any(x in line for x in ["Cp Coefficients", "Strip"])):
                break
            if data_start:
                try:
                    values = list(map(float, line.split()[:8]))
                    if len(values) == 8:
                        data.append(values)
                except ValueError:
                    continue
    return np.array(data[int(len(data)/2)::])

# Load AVL Data
xflr_file_0 = "WP4.1/XFLR0.txt"
xflr_data_0 = read_xflr_data(xflr_file_0)
xflr_file_10 = "WP4.1/XFLR10.txt"
xflr_data_10 = read_xflr_data(xflr_file_10)

# Extract Data
spanwise_positions = xflr_data_0[:, 0]
chords = xflr_data_0[:, 1]

Cls_0 = xflr_data_0[:, 3]
Cd_induced_0 = xflr_data_0[:, 5]

Cls_10 = xflr_data_10[:, 3]
Cd_induced_10 = xflr_data_10[:, 5]


Cm_0 = xflr_data_0[:, 7]
Cm_10 = xflr_data_10[:, 7]


#Cls_0= sp.interpolate.interp1d(spanwise_positions,Cls_0,kind='cubic',fill_value="extrapolate")
#Cls_10 = sp.interpolate.interp1d(spanwise_positions,Cls_10,kind='cubic',fill_value="extrapolate")



 # Calculate Distributed Loads
#L_dist = 0.5 * rho * velocity**2 * Cls_0_interpolated(spanwise_positions) * chords
#D_dist = 0.5 * rho * velocity**2 * Cd_induced_0 * chords
#N_dist = np.cos(alpha_a) * L_dist + np.sin(alpha_a) * D_dist


# # Engine Properties
engine_position = 3.9 # [m]
engine_weight = 2858 * 9.80665 #[N]
engine_torque = 240000  #[Nm]

# Load Factors
load_factor_positive = 2.5
load_factor_negative = -1.5

#distributed_load_positive = N_dist * load_factor_positive
#distributed_load_negative = N_dist * load_factor_negative


def coefficients(Cls0, Cls10, CLd, Cm_0, Cm_10 ):

    CLds = Cls0 + (CLd - Cls0)/(Cls10- Cls0) * ( Cls0- Cls10)
    alpha = (CLd - Cls0)/(Cls10-Cls0) * 10
    CD= CLds**2 / (np.pi * 8.05 * 0.891)
    CM= Cm_0 + (Cm_10- Cm_0)* alpha


    CN = CLds * np.cos(alpha* np.pi/180) +CD * np.sin(alpha* np.pi/180)
    CT = CLds * np.sin(alpha* np.pi/180) +CD * np.cos(alpha* np.pi/180)
    return(CN,CM)
CN, CM =coefficients(Cls_0, Cls_10, 0.5, Cm_0, Cm_10)

def dimensionalize(CN,CT,chords):
    rho = 1.225
    v = 225  #[m/s]
    N= CN * 0.5* rho* v**2 * chords
    T = CT* 0.5 * rho * v ** 2 * chords
    M = CM * 0.5* rho * v**2 * chords**2
    return(N,M)

N,M = dimensionalize(CN,10,chords)


def distributed_shear_force(normal, positions, point_loads):
    """
    Compute distributed shear force over a span of positions.
    
    Parameters:
        normal: Array or function of the distributed normal force
        positions: Array of spanwise positions
        lift: Total lift force to integrate up to
        point_loads: List of tuples (force, position) for discrete loads
    
    Returns:
        Shear force distribution as an array corresponding to `positions`.
    """
    shear_force = []
    for i, x in enumerate(positions):
        # Integral for distributed normal force from the current position
        integral, _ = quad(normal, x, max(spanwise_positions))

        # Contribution from point loads
        point_load_contribution = engine_weight
        if point_loads:
            for P, z_p in point_loads:
                if z_p >= x:  # Only include point loads outboard of x``
                    point_load_contribution += P

        # Shear force at current position
        S_x = -integral - point_load_contribution
        shear_force.append(S_x)
    return np.array(shear_force)

plt.plot(spanwise_positions, CN)
plt.show

#
# # Functions
# def interpolate_distributed_load(x, spanwise_positions, distributed_load):
#     return np.interp(x, spanwise_positions, distributed_load)
#
# def compute_shear_force(x_eval, spanwise_positions, distributed_load, point_load_position, point_load):
#     def distributed_load_function(x):
#         return interpolate_distributed_load(x, spanwise_positions, distributed_load)
#
#     integral_w, _ = quad(distributed_load_function, x_eval, wing_span)
#
#     S_eval = -integral_w
#     if x_eval <= point_load_position:
#         S_eval += point_load
#
#     return S_eval
#
# def bending_moment(x_eval, spanwise_positions, shear_force_function):
#     def shear_integral(x):
#         return shear_force_function(x)
#
#     integral_s, _ = quad(shear_integral, x_eval, wing_span)
#     return -integral_s
#
# def torque_distribution(x_eval, spanwise_positions, distributed_load, torque_position, torque_value):
#     def torque_integral(x):
#         interpolated_load = interpolate_distributed_load(x, spanwise_positions, distributed_load)
#         arm_length = 0.05  # DISTANCE FROM MIDDLE OF ENGINE TO WING LINE
#         return interpolated_load * arm_length
#
#     integral_torque, _ = quad(torque_integral, x_eval, wing_span)
#
#     T_eval = integral_torque
#     if x_eval <= torque_position:
#         T_eval += torque_value
#
#     return T_eval
#
# # Compute Results
# spanwise_positions2 = np.linspace(0, np.max(spanwise_positions), 1000)
#
# shear_force_positive = np.array(
#     [compute_shear_force(x, spanwise_positions, distributed_load_positive, engine_position, engine_weight)
#      for x in spanwise_positions2]
# )
# shear_force_negative = np.array(
#     [compute_shear_force(x, spanwise_positions, distributed_load_negative, engine_position, engine_weight)
#      for x in spanwise_positions2]
# )
#
# bending_moment_positive = np.array(
#     [bending_moment(x, spanwise_positions, lambda x: interpolate_distributed_load(x, spanwise_positions2, shear_force_positive))
#      for x in spanwise_positions2]
# )
# bending_moment_negative = np.array(
#     [bending_moment(x, spanwise_positions, lambda x: interpolate_distributed_load(x, spanwise_positions2, shear_force_negative))
#      for x in spanwise_positions2]
# )
#
# torque_positive = np.array(
#     [torque_distribution(x, spanwise_positions, distributed_load_positive, engine_position, engine_torque) for x in spanwise_positions2]
# )
# torque_negative = np.array(
#     [torque_distribution(x, spanwise_positions, distributed_load_negative, engine_position, engine_torque) for x in spanwise_positions2]
# )
#
# # Plot Results
# fig, axs = plt.subplots(3, 2, figsize=(15, 12))
#
# # Shear Force
# axs[0, 0].plot(spanwise_positions2, shear_force_positive, label="Shear Force (+)", color='blue')
# axs[0, 1].plot(spanwise_positions2, shear_force_negative, label="Shear Force (-)", color='red')
# axs[0, 0].set_title("Positive Load Factor - Shear Force")
# axs[0, 1].set_title("Negative Load Factor - Shear Force")
# axs[0, 0].set_ylabel("Shear Force (N)")
# axs[0, 1].set_ylabel("Shear Force (N)")
#
# # Bending Moment
# axs[1, 0].plot(spanwise_positions2, bending_moment_positive, label="Bending Moment (+)", color='blue')
# axs[1, 1].plot(spanwise_positions2, bending_moment_negative, label="Bending Moment (-)", color='red')
# axs[1, 0].set_title("Positive Load Factor - Bending Moment")
# axs[1, 1].set_title("Negative Load Factor - Bending Moment")
# axs[1, 0].set_ylabel("Bending Moment (Nm)")
# axs[1, 1].set_ylabel("Bending Moment (Nm)")
#
# # Torque
# axs[2, 0].plot(spanwise_positions2, torque_positive, label="Torque (+)", color='blue')
# axs[2, 1].plot(spanwise_positions2, torque_negative, label="Torque (-)", color='red')
# axs[2, 0].set_title("Positive Load Factor - Torque")
# axs[2, 1].set_title("Negative Load Factor - Torque")
# axs[2, 0].set_ylabel("Torque (Nm)")
# axs[2, 1].set_ylabel("Torque (Nm)")
#
# for ax in axs.flat:
#     ax.set_xlabel("Spanwise Position (m)")
#     ax.legend()
#     ax.grid()
#
# plt.tight_layout()
# plt.show()
#
# print(torque_positive.shape)


