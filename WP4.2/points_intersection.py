import numpy as np

def find_intersections(airfoil_points, x_lines):
    # print(airfoil_points)
    intersections = []
    airfoil_points = np.sort(airfoil_points, axis=0)
    airfoil_points = np.delete(airfoil_points,0,0)
    print(airfoil_points)
    listoflengths = np.array([[],[]])
    firstcolumn = airfoil_points[::2,0]
    secondcolumn = airfoil_points[::2,1] - airfoil_points[1::2,1]
    # for i in range(len(airfoil_points)):
    #     listoflengths += [airfoil_points[i+1][0], airfoil_points[i+1][1] - airfoil_points[i][1]]
    print("After: ", listoflengths)
    x_coords, y_coords = zip(*airfoil_points)
    return intersections

def run(x_lines):
    file_path = "WP4.2/fx60126.dat"
    airfoil_points = []
    airfoil_points = np.loadtxt(file_path)
    intersections = find_intersections(airfoil_points, x_lines)
    print(intersections)
    return intersections
        

