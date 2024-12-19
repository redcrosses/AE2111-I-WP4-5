def main5(I_xx, trapezoids, span_and_chord, loads, spanwise_position, max_stress=450e6, max_margin=10):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter1d

    def Mx(y, load_case):
        # Interpolate the moment along the spanwise position for the given load case
        moment = np.interp(y, spanwise_position, loads[load_case][1], left=0, right=0)
        
        # If the moment is exactly 0, replace it with the closest non-zero value
        if moment == 0:
            non_zero_values = np.array(loads[load_case][1])[np.array(loads[load_case][1]) != 0]
            closest_index = np.argmin(np.abs(np.array(spanwise_position) - y))
            moment = non_zero_values[np.argmin(np.abs(non_zero_values - loads[load_case][1][closest_index]))]
        
        # Ensure the moment meets the threshold requirement
        return moment

    # Initialize a figure for side-by-side subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    for load_case in range(2):  # Iterate through the two load cases
        margin_of_safety_list = []
        moment_list = []
        M_max = 0
        Failed = False

        for i in range(span_and_chord.shape[0]):
            I = I_xx[i]
            M = Mx(span_and_chord[i, 1], load_case)  # Moment at current spanwise position for load case
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

        # Plot the smoothed margin of safety for the current load case
        axes[load_case].plot(span_and_chord[:, 1], smoothed_margin, label="Smoothed Margin of Safety", color="blue", linewidth=2)
        axes[load_case].axhline(1, color="red", linestyle="--", label="Critical Safety Threshold")
        axes[load_case].set_xlabel("Spanwise Position [m]", fontsize=12)
        axes[load_case].set_ylabel("Margin of Safety [-]", fontsize=12)
        axes[load_case].set_title(f"Load Case {load_case + 1}", fontsize=14)
        axes[load_case].legend()
        axes[load_case].grid(True)

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

    return margin_of_safety_list.tolist(), M_max


if __name__ == "__main__":
    # Example inputs (replace with real data for testing)
    pass
