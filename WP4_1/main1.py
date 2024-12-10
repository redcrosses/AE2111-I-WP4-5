def main1(load_factor_1: float, load_factor_2: float,):
        
    import numpy as np
    import pandas as pd
    from scipy import interpolate, integrate
    import matplotlib.pyplot as plt
    from alive_progress import alive_bar

    # Constants
    rho = 0.412  # Air density [kg/m^3]
    V = 254.6  # Freestream velocity [m/s]
    q = 0.5 * rho * V**2  # Dynamic pressure [N/m^2]
    CL0 = 0.370799  # Lift coefficient at alpha = 0° (from file)
    CL10 = 1.171363  # Lift coefficient at alpha = 10° (from file)
    #D_x = 0.47  # Distance from shear center [m]
    w= 110389.610390
    s=134.81210

    # Function to process aerodynamic data
    def reprocess_aerodynamic_data(file_path):
        with open(file_path, 'r', encoding='latin1') as file:
            lines = file.readlines()

        # Locate the "Main Wing" section
        start_idx = None
        for i, line in enumerate(lines):
            if "Main Wing" in line:
                start_idx = i + 2
                break

        if start_idx is None:
            raise ValueError("Main Wing section not found in the file.")

        # Extract aerodynamic coefficients
        data = []
        for line in lines[start_idx:]:
            if line.strip() == "" or not line.lstrip().startswith(("-", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")):
                break
            split_line = line.split()
            try:
                y_span, chord, Cl, Cd, Cm = map(float, [split_line[0], split_line[1], split_line[3], split_line[5], split_line[6]])
                data.append([y_span, chord, Cl, Cd, Cm])
            except (ValueError, IndexError):
                continue

        return pd.DataFrame(data, columns=["y_span", "chord", "Cl", "Cd", "Cm"])

    # Function to compute angle of attack
    def compute_alpha(CL_d):
        return ((CL_d - CL0) / (CL10 - CL0)) * 10  # In degrees

    # Function to compute lift coefficient distribution
    def compute_cl_distribution(y_span, Cl_interp_a0, Cl_interp_a10, CL_d):
        Cl0_y = Cl_interp_a0(y_span)
        Cl10_y = Cl_interp_a10(y_span)
        return Cl0_y + ((CL_d - CL0) / (CL10 - CL0)) * (Cl10_y - Cl0_y)

    def compute_cd_distribution(y_span, Cd_interp_a0, Cd_interp_a10, CL_d):
        Cd0_y = Cd_interp_a0(y_span)
        Cd10_y = Cd_interp_a10(y_span)
        return Cd0_y + ((CL_d**2 - CL0**2) / (CL10**2 - CL0**2)) * (Cd10_y - Cd0_y)

    # Function to compute moment coefficient distribution
    def compute_cm_distribution(y_span, Cm_interp_a0, Cm_interp_a10, CL_d):
        Cm0_y = Cm_interp_a0(y_span)
        Cm10_y = Cm_interp_a10(y_span)
        return Cm0_y + ((CL_d - CL0) / (CL10 - CL0)) * (Cm10_y - Cm0_y)

    # Function to compute dimensional forces
    def compute_dimensional_forces(y_span, chord_interp, Cl_interp, Cd_interp, Cm_interp, q):
        chord = chord_interp(y_span)
        Cd = Cd_interp(y_span) + 0.0240256  # Corrected drag coefficient
        Cm = Cm_interp(y_span)

        L_prime = Cl_interp * q * chord  # Lift per unit span [N/m]
        D_prime = Cd * q * chord  # Drag per unit span [N/m]
        M_prime = Cm * q * chord**2  # Moment per unit span [Nm/m]

        return {"L_prime": L_prime, "D_prime": D_prime, "M_prime": M_prime}

    # Function to compute normal force distribution
    def compute_normal_force_distribution(L_prime, D_prime, alpha_d):
        alpha_rad = np.radians(alpha_d)  # Convert alpha to radians
        return np.cos(alpha_rad) * L_prime + np.sin(alpha_rad) * D_prime
    # Function to compute shear force distribution
    def compute_shear_force(y_span, N_prime, point_load=None, point_load_position=None):
        shear_force = []
        for i, y in enumerate(y_span):
            integral, _ = integrate.quad(lambda yp: np.interp(yp, y_span, N_prime), y, y_span[-1])
            S = integral  # Convention: upward forces are positive
            if point_load is not None and point_load_position is not None:
                if y <= point_load_position:
                    S += point_load
            shear_force.append(S)
        return np.array(shear_force)

    # Function to compute bending moment distribution
    def compute_bending_moment(y_span, shear_force, point_moment=None, point_moment_position=None):
        bending_moment = []
        for i, y in enumerate(y_span):
            integral, _ = integrate.quad(lambda yp: np.interp(yp, y_span, shear_force), y, y_span[-1])
            M = -integral  # Convention: clockwise moments are negative
            if point_moment is not None and point_moment_position is not None:
                if y <= point_moment_position:
                    M -= point_moment
            bending_moment.append(M)
        return np.array(bending_moment)


    def compute_dx_interpolator(y_span, chord_interp):

        # Compute D_x at the given y_span locations
        D_x_values = 0.47 - chord_interp(y_span) / 4

        # Create and return an interpolator for D_x
        return interpolate.interp1d(y_span, D_x_values, kind='cubic', fill_value="extrapolate")


    def compute_torque_distribution(y_span, N_prime, M_prime, D_x, point_load=None, point_load_position=None):
        torque = []
        for i, y in enumerate(y_span):
            q_torque = N_prime * D_x[i]+ M_prime  # Distributed torque per unit span
            integral, _ = integrate.quad(lambda yp: np.interp(yp, y_span, q_torque), y, y_span[-1])
            T = integral  # Integrate distributed torque
            if point_load is not None and point_load_position is not None:
                if y <= point_load_position:
                    T += point_load  # Use the specific D_x value at the current index
            torque.append(T)
        return np.array(torque)


    # Function to compute torque distribution
    #def compute_torque_distribution(y_span, N_prime, M_prime, D_x, point_load=None, point_load_position=None):
    #    torque = []
    #    for i, y in enumerate(y_span):
    #        q_torque = N_prime * D_x + M_prime  # Distributed torque per unit span
    #        integral, _ = integrate.quad(lambda yp: np.interp(yp, y_span, q_torque), y, y_span[-1])
    #        T = integral  # Integrate distributed torque
    #        if point_load is not None and point_load_position is not None:
    #            if y <= point_load_position:
    #                T += point_load * D_x
    #        torque.append(T)
    #    return np.array(torque)

    # Load aerodynamic data
    file_path_a0 = "WP4_1/XFLR0.txt"
    file_path_a10 = "WP4_1/XFLR10"

    #file_path_a0 = "MainWing_a=0.00_v=10.00ms.txt"
    #file_path_a10 = "MainWing_a=10.00_v=10.00ms.txt"

    df_a0 = reprocess_aerodynamic_data(file_path_a0)
    df_a10 = reprocess_aerodynamic_data(file_path_a10)

    # Filter for positive y_span
    df_a0_positive = df_a0[df_a0["y_span"] > 0]
    df_a10_positive = df_a10[df_a10["y_span"] > 0]

    # Interpolation functions
    Cl_interp_a0 = interpolate.interp1d(df_a0_positive["y_span"], df_a0_positive["Cl"], kind='cubic', fill_value="extrapolate")
    Cd_interp_a0 = interpolate.interp1d(df_a0_positive["y_span"], df_a0_positive["Cd"], kind='cubic', fill_value="extrapolate")
    Cm_interp_a0 = interpolate.interp1d(df_a0_positive["y_span"], df_a0_positive["Cm"], kind='cubic', fill_value="extrapolate")
    chord_interp_a0 = interpolate.interp1d(df_a0_positive["y_span"], df_a0_positive["chord"], kind='cubic', fill_value="extrapolate")
    Cl_interp_a10 = interpolate.interp1d(df_a10_positive["y_span"], df_a10_positive["Cl"], kind='cubic', fill_value="extrapolate")
    Cm_interp_a10 = interpolate.interp1d(df_a10_positive["y_span"], df_a10_positive["Cm"], kind='cubic', fill_value="extrapolate")
    Cd_interp_a10 = interpolate.interp1d(df_a10_positive["y_span"], df_a10_positive["Cd"], kind='cubic', fill_value="extrapolate")

    # Evaluation points
    y_span_eval = np.linspace(df_a0_positive["y_span"].min(), df_a0_positive["y_span"].max(), 1000)

    # Load cases
    load_cases = {"Positive Load Factor (n=2.5)": load_factor_1, "Negative Load Factor (n=-1.0)": load_factor_2}

    # Dictionary to store results
    results = {}

    # Iterate through load cases
    with alive_bar(6, title= "\033[96m {} \033[00m".format("WP4.1:"), bar='smooth', spinner='classic') as bar:
        for label, load_factor in load_cases.items():
            print("New Load Case")
            # Compute the desired lift coefficient (CL_d) based on the load factor
            CL_d = CL0*load_factor #(load_factor * w) / (q*s)
            alpha_d = compute_alpha(CL_d)
            # print(CL_d, CL0*load_factor)
            D_x_interpolator = compute_dx_interpolator(y_span_eval, chord_interp_a0)
            D_x = D_x_interpolator(y_span_eval)  # Evaluate D_x at any spanwise location

            # Compute aerodynamic distributions
            Cl_d_y = compute_cl_distribution(y_span_eval, Cl_interp_a0, Cl_interp_a10, CL_d)
            Cm_d_y = compute_cm_distribution(y_span_eval, Cm_interp_a0, Cm_interp_a10, CL_d)
            Cd_d_y = compute_cm_distribution(y_span_eval, Cd_interp_a0, Cd_interp_a10, CL_d)
            dimensional_forces = compute_dimensional_forces(
                y_span_eval, chord_interp_a0, Cl_d_y, lambda y: Cd_d_y, lambda y: Cm_d_y, q
            )
            N_prime = compute_normal_force_distribution(dimensional_forces["L_prime"], dimensional_forces["D_prime"], alpha_d)

            # Compute shear force, bending moment, and torque distributions
            print("Shear force distribution...")
            shear_force_distribution = compute_shear_force(y_span_eval, N_prime, point_load=2858 * 9.80665, point_load_position=3.9)
            bar()
            print("Bending moment distribution...")
            bending_moment_distribution = compute_bending_moment(y_span_eval, shear_force_distribution)
            bar()
            print("Torque distribution...")
            torque_distribution = compute_torque_distribution(
                y_span_eval, N_prime, dimensional_forces["M_prime"], D_x, point_load=240000, point_load_position=3.9
            )
            bar()
            
            # Store results for plotting
            results[label] = {
                "shear_force": shear_force_distribution,
                "bending_moment": bending_moment_distribution,
                "torque": torque_distribution
            }

    # Create a single figure for all 6 plots
    fig, axs = plt.subplots(3, 2, figsize=(15, 18))

    # Loop through results and plot on separate axes
    for i, (label, res) in enumerate(results.items()):
        # Shear Force
        axs[0, i].plot(y_span_eval, res["shear_force"], label=f"Shear Force - {label}")
        axs[0, i].axvline(x=3.9, color='r', linestyle='--', label="Point Load Position")
        axs[0, i].set_title(f"Shear Force - {label}")
        axs[0, i].set_xlabel("Spanwise Position (y) [m]")
        axs[0, i].set_ylabel("Shear Force [N]")
        axs[0, i].legend()
        axs[0, i].grid()

        # Bending Moment
        axs[1, i].plot(y_span_eval, res["bending_moment"], label=f"Bending Moment - {label}")
        axs[1, i].axvline(x=3.9, color='r', linestyle='--', label="Point Load Position")
        axs[1, i].set_title(f"Bending Moment - {label}")
        axs[1, i].set_xlabel("Spanwise Position (y) [m]")
        axs[1, i].set_ylabel("Bending Moment [Nm]")
        axs[1, i].legend()
        axs[1, i].grid()

        # Torque
        axs[2, i].plot(y_span_eval, res["torque"], label=f"Torque - {label}")
        axs[2, i].axvline(x=3.9, color='r', linestyle='--', label="Point Load Position")
        axs[2, i].set_title(f"Torque - {label}")
        axs[2, i].set_xlabel("Spanwise Position (y) [m]")
        axs[2, i].set_ylabel("Torque [Nm]")
        axs[2, i].legend()
        axs[2, i].grid()

    # Adjust layout
    fig.tight_layout()
    plt.show(block = False)
    
    results_pos = [
        list(results["Positive Load Factor (n=2.5)"]["shear_force"]),
        list(results["Positive Load Factor (n=2.5)"]["bending_moment"]),
        list(results["Positive Load Factor (n=2.5)"]["torque"])
    ]

    results_neg = [
        list(results["Negative Load Factor (n=-1.0)"]["shear_force"]),
        list(results["Negative Load Factor (n=-1.0)"]["bending_moment"]),
        list(results["Negative Load Factor (n=-1.0)"]["torque"])
    ]

    return results_pos, results_neg, y_span_eval

if __name__ == "__main__":
    results1, results2, results3 = main1(2.0, -1.5)
    print(results1[1])
    # print("Res pos:", results[0])
    # print("Res neg:", results[1])


    
