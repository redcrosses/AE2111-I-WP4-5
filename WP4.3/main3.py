import math
import matplotlib.pyplot as plt
import numpy as np
import intersect

# Constants and given values
a = 301.83  # Speed of sound at cruise altitude (m/s)
W = 170204.185  # Aircraft weight (lb)
Cl_clean = 1.44  # Max lift coefficient (clean configuration)
Cl_extended = 2.73  # Max lift coefficient (flaps extended)
rho = 0.441653  # Air density at cruise altitude (kg/m³) 1.55
S = 93.67  # Wing area (m²)
h = 9449.8  # Cruise altitude (m)

# Convert weight to Newtons for SI units
W_newtons = W * 4.44822  # 1 lb = 4.44822 N

# Stall speeds
Vs = math.sqrt((2 * W_newtons) / (rho * S * Cl_extended))  # Stall speed with flaps extended
Vs0 = math.sqrt((2 * W_newtons) / (rho * S * Cl_clean))    # Stall speed clean configuration

# Key speeds
Vc = 256.55  # Cruise speed (m/s)
Vd = 320.69  # Dive speed (m/s)

# Load factor limits
n_min = -1  # Minimum load factor
n_max_fd = 2  # Max load factor (flaps down)
n_max = 2.1 + (24000 / (W + 10000))  # Max load factor in clean configuration
if n_max<2.5:
	n_max=2.5
elif n_max>3.8:
	n_max=3.8


# Define velocity range for plotting
V = np.linspace(0, Vd, 500)  # Speeds from 0 to Vd
n_positive = (V / Vs0) ** 2  # Positive load factor curve (clean configuration)
n_positive[n_positive > n_max] = n_max  # Limit the curve to n_max

# Flaps-down region
n_positive_fd= (V/Vs)**2 
n_positive_fd[n_positive_fd > n_max_fd] = n_max_fd

# Negative load factor curve
n_negative = -n_positive
n_negative[n_negative < n_min] = n_min  # Limit the negative curve to n_min

#linear part from n_min to zero
# Linear segment from n_min to 0 between V_c and V_d
V_linear = np.linspace(Vc, Vd, 100)
n_linear = np.linspace(n_min, 0, 100)  # Linear interpolation from n_min to 0
# Plotting
plt.figure(figsize=(10, 6))

# Positive load factor curve
plt.plot(V, n_positive, label="Positive Load Factor", color="blue")
plt.plot(V,n_positive_fd, label="Positive Load Factor w. extended flaps", color="yellow")
# Negative load factor curve
plt.plot(V, n_negative, label="Negative Load Factor", color="green")


# Vertical lines for key speeds

plt.axvline(Vd, color="red", label="Dive Speed $V_d$")

# Horizontal lines for max/min load factors


plt.plot(V_linear, n_linear, color="red")

# Labels and legend
plt.title("V-n Diagram for Maneuvering Loads")
plt.xlabel("Velocity (m/s)")
plt.ylabel("Load Factor (n)")
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.ylim(n_min - 0.5, n_max + 0.5)
plt.xlim(0, Vd + 10)
plt.show()


plt.legend(loc="center left", bbox_to_anchor =(1.05,0.5))