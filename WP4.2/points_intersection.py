import numpy as np
def find_intersections(airfoil_points, y_vals):
    original = airfoil_points
    airfoil_points = np.array(airfoil_points)
    intersections = []
    airfoil_points = airfoil_points[np.argsort(airfoil_points[:,0])]
    airfoil_points = np.delete(airfoil_points,0,0)
    pos = airfoil_points[::2,0]
    lengths = np.abs(airfoil_points[1::2,1] - airfoil_points[::2,1])
    first = True
    for i in range(len(pos)):
        #front spar
        if first:
            try:
                if lengths[i] <= y_vals[0] <= lengths[i+1]: #trip on ascent
                    frontspar11, frontspar21 = zip(*original[np.where(original == pos[i]),:][0].transpose())
                    frontspar12, frontspar22 = zip(*original[np.where(original == pos[i+1]),:][0].transpose())                
                    frontsparpos = np.interp(y_vals[0], [lengths[i], lengths[i+1]], [pos[i], pos[i+1]])
                    frontspary1 = np.interp(frontsparpos, [frontspar11[0], frontspar12[0]], [frontspar11[1], frontspar12[1]])
                    frontspary2 = np.interp(frontsparpos, [frontspar21[0], frontspar22[0]], [frontspar21[1], frontspar22[1]])
                    # print("sparpos: ", frontsparpos, y_vals[0])
                    frontspar1 = np.array([frontsparpos, frontspary1])
                    frontspar2 = np.array([frontsparpos, frontspary2])

                    # print("Front: ", frontspar1, frontspar2)
                    first = False
            except:
                print("\033[91mFront spar length too large\033[00m")
        #rear spar
        else:
            try:
                if lengths[i+1] <= y_vals[1] <= lengths[i]: #trip on descent
                    rearspar11, rearspar21 = zip(*original[np.where(original == pos[i]),:][0].transpose())
                    rearspar12, rearspar22 = zip(*original[np.where(original == pos[i+1]),:][0].transpose())                
                    rearsparpos = np.interp(y_vals[1], [lengths[i+1], lengths[i]], [pos[i+1], pos[i]])
                    rearspary1 = np.interp(rearsparpos, [rearspar11[0], rearspar12[0]], [rearspar11[1], rearspar12[1]])
                    rearspary2 = np.interp(rearsparpos, [rearspar21[0], rearspar22[0]], [rearspar21[1], rearspar22[1]])
                    rearspar1 = np.array([rearsparpos, rearspary1])
                    rearspar2 = np.array([rearsparpos, rearspary2])

                    # print("Rear: ", rearspar1, rearspar2)
                    break
            except:
                print("\033[91mRear spar length too large\033[00m")
    intersections = np.array([frontspar1, frontspar2, rearspar2, rearspar1])
    # print("Trapezoid points:", intersections)
    return intersections
def run(y_vals):
    file_path = "WP4.2/fx60126.dat"
    airfoil_points = []
    airfoil_points = np.loadtxt(file_path)
    intersections = find_intersections(airfoil_points, y_vals)
    # print(intersections)
    return intersections
        

