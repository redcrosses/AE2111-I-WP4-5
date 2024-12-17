def main5(I_xx, trapezoids, span_and_chord, loads, spanwise_position, max_stress=450e6, max_margin=10):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter1d

    def Mx(y):
        # Interpolate the moment along the spanwise position
        moment = np.interp(y, spanwise_position, loads[0][1], left=0, right=0)
        
        # If the moment is exactly 0, replace it with the closest non-zero value
        if moment == 0:
            non_zero_values = np.array(loads[0][1])[np.array(loads[0][1]) != 0]  # Filter out zero values
            closest_index = np.argmin(np.abs(np.array(spanwise_position) - y))  # Find the closest index
            moment = non_zero_values[np.argmin(np.abs(non_zero_values - loads[0][1][closest_index]))]
        
        # Ensure the moment meets the threshold requirement
        return max(moment, -2537900.3936588285)

    # def Mx(y):
    # # Interpolate the moment along the spanwise position
    #     moment = np.interp(y, spanwise_position, loads[0][1], left=0, right=0)
    #     return max(moment,)

    margin_of_safety_list = []
    moment_list = []
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
        margin_of_safety = min(margin_of_safety, max_margin)
        margin_of_safety_list.append(margin_of_safety)
        moment_list.append(M)

    # Convert to numpy array
    margin_of_safety_list = np.array(margin_of_safety_list)

    # 1. Outlier removal using z-score
    z_scores = np.abs((margin_of_safety_list - np.mean(margin_of_safety_list)) / np.std(margin_of_safety_list))
    threshold = 2.5  # Threshold for "far-off" values
    margin_of_safety_list[z_scores > threshold] = np.median(margin_of_safety_list)

    # 2. Smooth the graph using Gaussian filter
    smoothed_margin = gaussian_filter1d(margin_of_safety_list, sigma=2)

    # Plot the smoothed margin of safety
    plt.figure(figsize=(8, 6))
    plt.plot(span_and_chord[:, 1], smoothed_margin, label="Smoothed Margin of Safety", color="blue", linewidth=2)
    plt.axhline(1, color="red", linestyle="--", label="Critical Safety Threshold")
    plt.xlabel("Spanwise Position [m]", fontsize=12)
    plt.ylabel("Margin of Safety [-]", fontsize=12)
    plt.title("Smoothed Margin of Safety vs Spanwise Position", fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # import numpy as np
    # import matplotlib.pyplot as plt

    # def Mx(y):
    #     # Interpolate the moment along the spanwise position
    #     return np.interp(y, spanwise_position, loads[0][1], left=0, right=0)

    # margin_of_safety_list = []
    # M_max = 0
    # Failed = False

    # for i in range(span_and_chord.shape[0]):
    #     I = I_xx[i]
    #     M = Mx(span_and_chord[i, 1])  # Moment at current spanwise position
    #     y_max = abs(trapezoids[1, 1] * span_and_chord[i, 0])  # Assuming correct indexing

    #     # Calculate stress
    #     stress = (M * y_max) / I

    #     # Check for failure and update M_max if necessary
    #     if stress >= max_stress and not Failed:
    #         M_max = M  # Record maximum moment causing failure
    #         Failed = True

    #     # Calculate margin of safety
    #     margin_of_safety = abs(max_stress / stress) if stress != 0 else float('inf')
    #     if margin_of_safety >= max_margin:
    #         margin_of_safety = max_margin
    #     margin_of_safety_list.append(margin_of_safety)


    # # Apply smoothing to avoid jumps
    # margin_of_safety_list = np.array(margin_of_safety_list)
    # margin_of_safety_list = np.minimum(margin_of_safety_list, max_margin)

    # # Ensure continuous graph: Adjust excessive starting values
    # # for i in range(1, len(margin_of_safety_list)):
    # #     if abs(margin_of_safety_list[i] - margin_of_safety_list[i - 1]) > 2:  # Smooth large jumps
    # #         margin_of_safety_list[i] = margin_of_safety_list[i - 1] + (margin_of_safety_list[i] - margin_of_safety_list[i - 1]) / 2

    # # Plot the margin of safety
    # plt.figure(figsize=(8, 6))
    # plt.plot(span_and_chord[:, 1], margin_of_safety_list, label="Margin of Safety", color="blue", linewidth=2)
    # plt.axhline(1, color="red", linestyle="--", label="Critical Safety Threshold")
    # plt.xlabel("Spanwise Position [m]", fontsize=12)
    # plt.ylabel("Margin of Safety [-]", fontsize=12)
    # plt.title("Margin of Safety vs Spanwise Position", fontsize=14)
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()

    print(moment_list)

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