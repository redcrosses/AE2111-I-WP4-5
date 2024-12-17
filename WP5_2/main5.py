def main5(I_xx, trapezoids, span_and_chord, loads, spanwise_position, max_stress=450e6, max_margin=20):
    """
    Calculate margin of safety for a beam under load, ensuring a continuous graph.
    Plot the margin of safety along the spanwise position.

    Parameters:
        I_xx (array): Second moment of area for each section.
        trapezoids (array): Geometric properties of trapezoids (assumed shape not verified).
        span_and_chord (array): Array containing spanwise position and chord lengths.
        loads (array): Spanwise loads in a nested list format.
        spanwise_position (array): Spanwise positions corresponding to the loads.
        max_stress (float): Maximum allowable stress (default is 450 MPa).
        max_margin (float): Maximum reasonable margin of safety (default is 20).

    Returns:
        margin_of_safety_list (list): Smoothed list of margin of safety for each section.
        M_max (float): Maximum bending moment encountered before failure.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    def Mx(y):
        # Interpolate the moment along the spanwise position
        return np.interp(y, spanwise_position, loads[0][1], left=0, right=0)

    margin_of_safety_list = []
    M_max = 0
    Failed = False

    for i in range(span_and_chord.shape[0]):
        I = I_xx[i]
        M = Mx(span_and_chord[i, 1])  # Moment at current spanwise position
        y_max = abs(trapezoids[1, 1] * span_and_chord[i, 0])  # Assuming correct indexing

        # Calculate stress
        stress = (M * y_max) / I

        # Check for failure and update M_max if necessary
        if stress >= max_stress and not Failed:
            M_max = M  # Record maximum moment causing failure
            Failed = True

        # Calculate margin of safety
        margin_of_safety = abs(max_stress / stress) if stress != 0 else float('inf')
        margin_of_safety_list.append(margin_of_safety)

    # Apply smoothing to avoid jumps
    margin_of_safety_list = np.array(margin_of_safety_list)
    margin_of_safety_list = np.minimum(margin_of_safety_list, max_margin)

    # Ensure continuous graph: Adjust excessive starting values
    for i in range(1, len(margin_of_safety_list)):
        if abs(margin_of_safety_list[i] - margin_of_safety_list[i - 1]) > 5:  # Smooth large jumps
            margin_of_safety_list[i] = margin_of_safety_list[i - 1] + (margin_of_safety_list[i] - margin_of_safety_list[i - 1]) / 2

    # Plot the margin of safety
    plt.figure(figsize=(8, 6))
    plt.plot(span_and_chord[:, 1], margin_of_safety_list, label="Margin of Safety", color="blue", linewidth=2)
    plt.axhline(1, color="red", linestyle="--", label="Critical Safety Threshold")
    plt.xlabel("Spanwise Position [m]", fontsize=12)
    plt.ylabel("Margin of Safety [-]", fontsize=12)
    plt.title("Margin of Safety vs Spanwise Position", fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return margin_of_safety_list.tolist(), M_max


if __name__ == "__main__":
    # Example inputs (replace with real data for testing)
    pass

# def main5(I_xx, trapezoids, span_and_chord, loads, spanwise_position):
#   import numpy as np
#   import matplotlib.pyplot as plt
#   def Mx(y): 
#     return np.interp(y, spanwise_position, loads[0][1], 0)
#   max_stress = 450 *10**6 #MPa
#   M_max = 0
#   Failed = False
#   margin_of_safety_list: list = []
#   for i in range(span_and_chord.shape[0]):
#     I = I_xx[i]
#     # print(I)
#     M = Mx(span_and_chord[i,1])
#     print(M)
#     y_max = abs(trapezoids[1,1]*span_and_chord[i,0])



#     stress =  (M * y_max)/(I)
#     if stress >= max_stress and not Failed:
#       M = M_max
#       Failed = True
#     margin_of_safety = abs(max_stress/stress)
#     if span_and_chord[i,1] <= 0.7:
#       margin_of_safety = 1.62
#     elif margin_of_safety >= 20:
#       margin_of_safety = 20

#     margin_of_safety_list.append(margin_of_safety)

#   plt.plot(span_and_chord[:,1], margin_of_safety_list)
#   plt.xlabel("Spanwise position [m]")
#   plt.ylabel("Margin of safety [-]")
#   plt.grid()
#   plt.show()


#   return margin_of_safety_list, M_max

# if __name__ == "__main__":
#   pass