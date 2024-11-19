def centroid_of_quadrilateral(points):
    """
    Calculate the centroid of a quadrilateral given its vertices.
    
    :param points: List of tuples [(x1, y1), (x2, y2), (x3, y3), (x4, y4)] in counterclockwise order.
    :return: Centroid coordinates (Cx, Cy).
    """
    # Ensure the quadrilateral is closed (last point = first point)
    points.append(points[0])
    
    # Calculate area A
    A = 0
    for i in range(len(points) - 1):
        A += points[i][0] * points[i + 1][1] - points[i + 1][0] * points[i][1]
    A = A / 2.0
    
    # Calculate centroid (Cx, Cy)
    Cx, Cy = 0, 0
    for i in range(len(points) - 1):
        factor = (points[i][0] * points[i + 1][1] - points[i + 1][0] * points[i][1])
        Cx += (points[i][0] + points[i + 1][0]) * factor
        Cy += (points[i][1] + points[i + 1][1]) * factor
    
    Cx /= (6 * A)
    Cy /= (6 * A)
    
    return Cx, Cy

# Example usage
quadrilateral = [(0, 0), (4, 0), (4, 3), (0, 6)]  # Example coordinates
centroid = centroid_of_quadrilateral(quadrilateral)
print("Centroid:", centroid)
