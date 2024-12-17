def stringer_sizing(length_1,length_2, thickness_1, thickness_2):
        # t1 and l1 are the dimensions of the rectangle parallel to the wingbox and t2 and l2 are the dimensions for the rectangle,
        #perperndicular to the wing box, and d1 is the distance between the x axis and the centroid of the parallel rectangle, 
        #higher order terms cancel out, MADE BY VICTOR BOSS
        Area_parallel = thickness_1 * length_1
        Area_perpendicular = thickness_2 * length_2 
        Total_area_stringer = Area_parallel + Area_perpendicular
        y_centroid_stringer = (Area_perpendicular*length_2/2)/Total_area_stringer
        x_centroid_stringer = (Area_parallel*length_1/2)/Total_area_stringer
        Ixx_stringer =  Area_parallel * y_centroid_stringer ** 2 + (1/12) * length_2 ** 3 * thickness_2 + Area_perpendicular * (length_2-y_centroid_stringer) ** 2
        Iyy_stringer =  Area_perpendicular * x_centroid_stringer ** 2 + (1/12) * length_1 ** 3 * thickness_1 + Area_parallel * (length_1-x_centroid_stringer) ** 2

        return Ixx_stringer, Iyy_stringer, Area_perpendicular, Area_parallel, Total_area_stringer
