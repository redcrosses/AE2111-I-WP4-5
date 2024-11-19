import numpy as np
import pandas as pd  # Ensure pandas is imported

def find_intersections(airfoil_points, x_lines):
    """
    Find the intersection points of the airfoil shape with the vertical bounding lines of the wingbox.
    
    Parameters:
    - airfoil_points: List of tuples (x, y) representing the airfoil shape.
    - x_lines: List of x-coordinates [x1, x2] representing the vertical bounding lines of the wingbox.
    
    Returns:
    - List of tuples (x, y) representing intersection points for each x_line.
    """
    intersections = []
    x_coords, y_coords = zip(*airfoil_points)
    
    for x_line in x_lines:
        # Find indices where the airfoil crosses the vertical line
        for i in range(len(x_coords) - 1):
            if (x_coords[i] <= x_line <= x_coords[i+1]) or (x_coords[i+1] <= x_line <= x_coords[i]):
                # Interpolate to find the y-coordinate
                x1, y1 = x_coords[i], y_coords[i]
                x2, y2 = x_coords[i+1], y_coords[i+1]
                y_intersection = y1 + (y2 - y1) * (x_line - x1) / (x2 - x1)
                intersections.append((x_line, y_intersection))
                break
    return intersections

# Read airfoil points from the given file path
file_path = "WP4.2/fx60126.dat"

# Initialize lists to store airfoil points
airfoil_points = []

try:
    # Open and read the file
    with open(file_path, "r") as f:
        lines = f.readlines()
    
    # Skip the first line (header)
    for line in lines[1:]:
        # Parse the x, y points from the file
        point = line.strip().split()
        if len(point) == 2:
            x, y = map(float, point)
            airfoil_points.append((x, y))
except FileNotFoundError:
    print(f"File not found: {file_path}")

# If points were read, calculate intersections
if airfoil_points:
    # Define the wingbox bounding lines
    x_lines = [0.1, 0.5]

    # Calculate intersections
    intersections = find_intersections(airfoil_points, x_lines)
    
    # Convert results to a DataFrame for visualization
    intersections_df = pd.DataFrame(intersections, columns=["x", "y"])
    print(intersections_df)  # Display the DataFrame instead of using undefined tools
else:
    print("No airfoil points were loaded from the file.")
