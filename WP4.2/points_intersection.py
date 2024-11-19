import numpy as np

def find_intersections(airfoil_points, x_lines):
    # print(airfoil_points)
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

def run(x_lines):
    file_path = "WP4.2/fx60126.dat"
    airfoil_points = []
    airfoil_points = np.loadtxt(file_path)
    intersections = find_intersections(airfoil_points[0:int(len(airfoil_points)/2)], x_lines)
    intersections2 = find_intersections(airfoil_points[int(len(airfoil_points)/2):], x_lines)
    intersections += intersections2[::-1]
    print(intersections)
    return intersections
        

