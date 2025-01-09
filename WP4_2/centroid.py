import numpy as np
def centroid_of_quadrilateral(numpy_list):
    points = [tuple(row) for row in numpy_list]
   
    # Ensure the polygon is closed by appending the first point to the end
    points.append(points[0])

    # Initialize variables for the centroid calculation
    signed_area = 0
    centroid_x = 0
    centroid_y = 0

    # Iterate over the edges of the polygon
    for i in range(len(points) - 1):
        x0, y0 = points[i]
        x1, y1 = points[i + 1]
        
        # Calculate the cross product (determinant) of the edge
        cross = x0 * y1 - x1 * y0
        signed_area += cross
        
        # Accumulate the centroid components
        centroid_x += (x0 + x1) * cross
        centroid_y += (y0 + y1) * cross

    # Finalize the calculations
    signed_area *= 0.5
    centroid_x /= (6 * signed_area)
    centroid_y /= (6 * signed_area)

    return centroid_x, centroid_y

if __name__ == "__main__":
    coords = np.array([[0, 1], [0, 0], [1, 0], [1, 1]])
    print(centroid_of_quadrilateral(coords))
