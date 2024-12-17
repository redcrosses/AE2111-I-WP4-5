import numpy as np

def critical_shear_stress(ks, t, b, E = 72.4*10**9, nu = 0.33, ):
    """
    Calculate the critical web buckling shear stress (τ_cr)

    Parameters:
        E: Young's modulus (Elastic modulus) in Pascals or consistent units
        nu: Poisson's ratio (dimensionless)
        ks: Shear buckling coefficient, depends on aspect ratio a/b (dimensionless)
        t: Thickness of the spar in meters
        b: Short side of the plate in meters

    Returns:
        Critical shear buckling stress τ_cr in Pa
    """

    # Calculate the critical shear buckling stress using NumPy operations
    tau_cr = (np.pi ** 2 * ks * E) / (12 * (1 - nu ** 2)) * (t / b) ** 2

    return tau_cr

# Example usage:
print(critical_shear_stress(5, 0.001, 0.5))
