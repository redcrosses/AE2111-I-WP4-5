import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Constants
wing_span = 30  # Total span of the wing (m)
alpha_a = 10 * np.pi / 180  # Angle of attack in radians
rho = 1.225  # Air density in kg/m^3 (sea level standard)
velocity = 10  # Freestream velocity in m/s

# Read AVL data from file
def read_avl_data():
    data_start = False
    data = []
    with open('AVL.txt', 'r') as f:
        for line in f:
            # Detect the start of numerical data
            if "y-span" in line:
                data_start = True
                continue
            if data_start:
                try:
                    values = list(map(float, line.split()[:5]))  # First 5 columns of interest
                    data.append(values)
                except ValueError:
                    break  # Stop reading at the end of data
    return np.array(data)

# Process AVL data
avl_file = "AVL.txt"  # Your file is named "AVL"
avl_data = read_avl_data()

# Extract spanwise position (y), chord (c), Cl, and Cd
spanwise_positions = avl_data[:, 0]  # y-span
chords = avl_data[:, 1]  # Chord
Cls = avl_data[:, 3]  # Cl
Cds = avl_data[:, 4]  # Cd

# Calculate distributed lift and drag per unit span
L_dist = 0.5 * rho * velocity ** 2 * Cls * chords  # Lift distribution (N/m)
D_dist = 0.5 * rho * velocity ** 2 * Cds * chords  # Drag distribution (N/m)

# Transform aerodynamic forces to body reference frame
N_dist = np.cos(alpha_a) * L_dist + np.sin(alpha_a) * D_dist  # Normal force distribution

# Set up engine properties
engine_position = 10  # Distance of engine from root (m)
engine_weight = 2000  # Engine weight (N)
engine_torque = 10000  # Engine torque (Nm)

# Calculate distributed load for positive and negative load factors
load_factor_positive = 2.5
load_factor_negative = -1.5
distributed_load_positive = N_dist * load_factor_positive
distributed_load_negative = N_dist * load_factor_negative

# Shear Force Calculation
def shear_force(x_eval, distributed_load, point_load_position, point_load):
    """Compute the shear force using numerical integration."""

    def integrand(x):
        idx = np.abs(spanwise_positions - x).argmin()
        return distributed_load[idx]

    S_eval = -integrate.quad(integrand, x_eval, wing_span)[0]
    if x_eval <= point_load_position:
        S_eval -= point_load
    return S_eval

shear_force_positive = np.array(
    [shear_force(x, distributed_load_positive, engine_position, engine_weight) for x in spanwise_positions])
shear_force_negative = np.array(
    [shear_force(x, distributed_load_negative, engine_position, engine_weight) for x in spanwise_positions])

# Bending Moment Calculation
def bending_moment(x_eval, shear_force_values):
    """Compute the bending moment using numerical integration."""

    def integrand(x):
        idx = np.abs(spanwise_positions - x).argmin()
        return shear_force_values[idx]

    M_eval = -integrate.quad(integrand, x_eval, wing_span)[0]
    return M_eval

bending_moment_positive = np.array([bending_moment(x, shear_force_positive) for x in spanwise_positions])
bending_moment_negative = np.array([bending_moment(x, shear_force_negative) for x in spanwise_positions])

# Torque Distribution
def torque_distribution(x_eval, distributed_load, torque_position, torque_value):
    """Compute torque distribution using numerical integration."""

    def integrand(x):
        idx = np.abs(spanwise_positions - x).argmin()
        return distributed_load[idx] * 0.05  # Assume constant moment arm of 0.05m

    T_eval = integrate.quad(integrand, x_eval, wing_span)[0]
    if x_eval <= torque_position:
        T_eval += torque_value
    return T_eval

torque_positive = np.array(
    [torque_distribution(x, distributed_load_positive, engine_position, engine_torque) for x in spanwise_positions])
torque_negative = np.array(
    [torque_distribution(x, distributed_load_negative, engine_position, engine_torque) for x in spanwise_positions])

# Plot Results
fig, axs = plt.subplots(3, 2, figsize=(15, 12))

# Shear Force Diagrams
axs[0, 0].plot(spanwise_positions, shear_force_positive, label="Shear Force (+)", color='blue')
axs[0, 1].plot(spanwise_positions, shear_force_negative, label="Shear Force (-)", color='red')
axs[0, 0].set_title("Positive Load Factor - Shear Force")
axs[0, 1].set_title("Negative Load Factor - Shear Force")
axs[0, 0].set_ylabel("Shear Force (N)")
axs[0, 1].set_ylabel("Shear Force (N)")

# Bending Moment Diagrams
axs[1, 0].plot(spanwise_positions, bending_moment_positive, label="Bending Moment (+)", color='blue')
axs[1, 1].plot(spanwise_positions, bending_moment_negative, label="Bending Moment (-)", color='red')
axs[1, 0].set_title("Positive Load Factor - Bending Moment")
axs[1, 1].set_title("Negative Load Factor - Bending Moment")
axs[1, 0].set_ylabel("Bending Moment (Nm)")
axs[1, 1].set_ylabel("Bending Moment (Nm)")

# Torque Diagrams
axs[2, 0].plot(spanwise_positions, torque_positive, label="Torque (+)", color='blue')
axs[2, 1].plot(spanwise_positions, torque_negative, label="Torque (-)", color='red')
axs[2, 0].set_title("Positive Load Factor - Torque")
axs[2, 1].set_title("Negative Load Factor - Torque")
axs[2, 0].set_ylabel("Torque (Nm)")
axs[2, 1].set_ylabel("Torque (Nm)")

# Common labels
for ax in axs.flat:
    ax.set_xlabel("Spanwise Position (m)")
    ax.legend()
    ax.grid()

plt.tight_layout()
plt.show()
