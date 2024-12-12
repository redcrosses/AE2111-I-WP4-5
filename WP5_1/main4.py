# t1 and l1 are the dimensions of the rectangle parallel to the wingbox and t2 and l2 are the dimensions for the rectangle
#perperndicular to the wing box
def stringer_sizing(length_1,length_2, thickness_1, thickness_2,):
    Ixx =  length_1 * thickness_1 * d1 ** 2 + (1/12) * length_2 ** 3 * thickness_2
    Area_parallel = thickness_1 * length_1
    Area_perpendicular = thickness_2 * length_2 
    Total_area_stringer = Area_parallel + Area_perpendicular
    return Ixx, Area_perpendicular, Area_parallel, Total_area_stringer