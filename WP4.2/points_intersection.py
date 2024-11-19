import numpy as np

def find_intersections(airfoil_points, y_vals):
    # print(airfoil_points)
    airfoil_points = np.array(airfoil_points)
    intersections = []
    airfoil_points = airfoil_points[np.argsort(airfoil_points[:,0])]
    airfoil_points = np.delete(airfoil_points,0,0)
    pos = airfoil_points[::2,0]
    lengths = np.abs(airfoil_points[1::2,1] - airfoil_points[::2,1])
    print(pos, lengths)
    # try:
        #front spar
    for i in range(len(pos)):
        print(lengths[i], y_vals[0]/2,lengths[i+1])
        if lengths[i] <= y_vals[0]/2 <= lengths[i+1]: #trip on ascent
            frontspar1 = airfoil_points[i-1]
            frontspar2 = airfoil_points[i]
            break

    #rear spar
    for i in range(len(pos)):
        if lengths[i+1] <= y_vals[1]/2 <= lengths[i]: #trip on descent
            rearspar1 = airfoil_points[i-1] 
            rearspar2 = airfoil_points[i]
            break
    intersections = [frontspar1, frontspar2, rearspar2,rearspar1]
    return intersections

def run(y_vals):
    file_path = "WP4.2/fx60126.dat"
    airfoil_points = []
    airfoil_points = np.loadtxt(file_path)
    intersections = find_intersections(airfoil_points, y_vals)
    print(intersections)
    return intersections
        

