def main1():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import quad

    # Constants
    wing_span = 27.47721  # Total span of the wing (m)
    alpha_a = 14 * np.pi / 180  # Angle of attack in radians
    rho = 1.225  # Air density in kg/m^3 (sea level standard)
    velocity = 0.85 * 343  # Freestream velocity in m/s

    # Read AVL Data
    def read_avl_data(avl_file):
        data_start = False
        data = []
        with open(avl_file, 'r') as f:
            for line in f:
                if "Main Wing" in line and not data_start:
                    data_start = True
                    continue
                if data_start and (not line.strip() or any(x in line for x in ["Cp Coefficients", "Strip"])):
                    break
                if data_start:
                    try:
                        values = list(map(float, line.split()[:6]))
                        if len(values) == 6:
                            data.append(values)
                    except ValueError:
                        continue
        return np.array(data[int(len(data)/2)::])

    # Load AVL Data
    avl_file = "WP4_1/AVL.txt"
    avl_data = read_avl_data(avl_file)

    # Extract Data
    spanwise_positions = avl_data[:, 0]
    chords = avl_data[:, 1]
    Cls = avl_data[:, 3]
    Cds = avl_data[:, 5]

    # Calculate Distributed Loads
    L_dist = 0.5 * rho * velocity**2 * Cls * chords
    D_dist = 0.5 * rho * velocity**2 * Cds * chords
    N_dist = np.cos(alpha_a) * L_dist + np.sin(alpha_a) * D_dist

    # Engine Properties
    engine_position = 3.9
    engine_weight = 56016.8
    engine_torque = 240000

    # Load Factors
    load_factor_positive = 2.5
    load_factor_negative = -1.5

    distributed_load_positive = N_dist * load_factor_positive
    distributed_load_negative = N_dist * load_factor_negative

    # Functions
    def interpolate_distributed_load(x, spanwise_positions, distributed_load):
        return np.interp(x, spanwise_positions, distributed_load)

    def compute_shear_force(x_eval, spanwise_positions, distributed_load, point_load_position, point_load):
        def distributed_load_function(x):
            return interpolate_distributed_load(x, spanwise_positions, distributed_load)

        integral_w, _ = quad(distributed_load_function, x_eval, wing_span)

        S_eval = -integral_w
        if x_eval <= point_load_position:
            S_eval += point_load

        return S_eval

    def bending_moment(x_eval, spanwise_positions, shear_force_function):
        def shear_integral(x):
            return shear_force_function(x)

        integral_s, _ = quad(shear_integral, x_eval, wing_span)
        return -integral_s

    def torque_distribution(x_eval, spanwise_positions, distributed_load, torque_position, torque_value):
        def torque_integral(x):
            interpolated_load = interpolate_distributed_load(x, spanwise_positions, distributed_load)
            arm_length = 0.05  # DISTANCE FROM MIDDLE OF ENGINE TO WING LINE
            return interpolated_load * arm_length

        integral_torque, _ = quad(torque_integral, x_eval, wing_span)

        T_eval = integral_torque
        if x_eval <= torque_position:
            T_eval += torque_value

        return T_eval

    # Compute Results
    spanwise_positions2 = np.linspace(0, np.max(spanwise_positions), 1000)

    shear_force_positive = np.array(
        [compute_shear_force(x, spanwise_positions, distributed_load_positive, engine_position, engine_weight)
        for x in spanwise_positions2]
    )
    shear_force_negative = np.array(
        [compute_shear_force(x, spanwise_positions, distributed_load_negative, engine_position, engine_weight)
        for x in spanwise_positions2]
    )

    bending_moment_positive = np.array(
        [bending_moment(x, spanwise_positions, lambda x: interpolate_distributed_load(x, spanwise_positions2, shear_force_positive))
        for x in spanwise_positions2]
    )
    bending_moment_negative = np.array(
        [bending_moment(x, spanwise_positions, lambda x: interpolate_distributed_load(x, spanwise_positions2, shear_force_negative))
        for x in spanwise_positions2]
    )

    torque_positive = np.array(
        [torque_distribution(x, spanwise_positions, distributed_load_positive, engine_position, engine_torque) for x in spanwise_positions2]
    )
    torque_negative = np.array(
        [torque_distribution(x, spanwise_positions, distributed_load_negative, engine_position, engine_torque) for x in spanwise_positions2]
    )

    # Plot Results
    fig, axs = plt.subplots(3, 2, figsize=(15, 12))

    # Shear Force
    axs[0, 0].plot(spanwise_positions2, shear_force_positive, label="Shear Force (+)", color='blue')
    axs[0, 1].plot(spanwise_positions2, shear_force_negative, label="Shear Force (-)", color='red')
    axs[0, 0].set_title("Positive Load Factor - Shear Force")
    axs[0, 1].set_title("Negative Load Factor - Shear Force")
    axs[0, 0].set_ylabel("Shear Force (N)")
    axs[0, 1].set_ylabel("Shear Force (N)")

    # Bending Moment
    axs[1, 0].plot(spanwise_positions2, bending_moment_positive, label="Bending Moment (+)", color='blue')
    axs[1, 1].plot(spanwise_positions2, bending_moment_negative, label="Bending Moment (-)", color='red')
    axs[1, 0].set_title("Positive Load Factor - Bending Moment")
    axs[1, 1].set_title("Negative Load Factor - Bending Moment")
    axs[1, 0].set_ylabel("Bending Moment (Nm)")
    axs[1, 1].set_ylabel("Bending Moment (Nm)")

    # Torque
    axs[2, 0].plot(spanwise_positions2, torque_positive, label="Torque (+)", color='blue')
    axs[2, 1].plot(spanwise_positions2, torque_negative, label="Torque (-)", color='red')
    axs[2, 0].set_title("Positive Load Factor - Torque")
    axs[2, 1].set_title("Negative Load Factor - Torque")
    axs[2, 0].set_ylabel("Torque (Nm)")
    axs[2, 1].set_ylabel("Torque (Nm)")

    for ax in axs.flat:
        ax.set_xlabel("Spanwise Position (m)")
        ax.legend()
        ax.grid()

    plt.tight_layout()
    plt.show()

    print(torque_positive.shape)


