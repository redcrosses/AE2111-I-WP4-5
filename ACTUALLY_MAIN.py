from WP4_1.main1 import main1
from WP4_2.main2 import main2
from WP4_3.main3 import main3
from WP5_1.main4 import main4
from WP5_2.main5 import main5
import sys
import matplotlib.pyplot as plt

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

print("\033[96m {} \033[00m".format("WP4.3:"))
print("Finding load factors")
n_positive, n_negative = main3()
print("Finding planform loading diagrams...")
loads_positive, loads_negative, spanwise_position = main1(n_positive, n_negative)

#creating the design and testing deflections
design1: object = main2(1, (loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.1056927, 0.07702, 0.01, 0.01, 2, 10, 0.04, 0.04, 0.001, 0.001) #loads, span positions, loading factors,  front spar length, rear spar length, horizontal spar thickness, vertical spar thickness, rib spacing [m], number of stringers #last four are dimensions of an L stringer; width, height, width thickness, height thickness [m]
design1.graph_displacements()
margin, max = main5(design1.moi_x_list, design1.trapezoid,  design1.chords_along_span, (loads_positive, loads_negative), spanwise_position) #testing tensile side of the design
main4(design1.moi_x_list, design1.trapezoid, design1.stringers, design1.chords_along_span, (loads_positive, loads_negative), spanwise_position, design1) #testing critical buckling case for the design
design1.graph_visualized_and_dash()

design2: object = main2(2, (loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.07702, 0.01, 0.01, 2, 10, 0.04, 0.04, 0.001, 0.001) #loads, span positions, loading factors,  front spar length, rear spar length, horizontal spar thickness, vertical spar thickness, rib spacing [m], number of stringers #last four are dimensions of an L stringer; width, height, width thickness, height thickness [m]
design2.graph_displacements()
margin, max = main5(design2.moi_x_list, design2.trapezoid,  design2.chords_along_span, (loads_positive, loads_negative), spanwise_position) #testing tensile side of the design
main4(design2.moi_x_list, design2.trapezoid, design2.stringers, design2.chords_along_span, (loads_positive, loads_negative), spanwise_position, design2) #testing critical buckling case for the design
design2.graph_visualized_and_dash()

design3: object = main2(3, (loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.06808, 0.01, 0.015, 3, 10, 0.04, 0.04, 0.001, 0.001) #loads, span positions, loading factors,  front spar length, rear spar length, horizontal spar thickness, vertical spar thickness, rib spacing [m], number of stringers #last four are dimensions of an L stringer; width, height, width thickness, height thickness [m]
design3.graph_displacements()
margin, max = main5(design3.moi_x_list, design3.trapezoid,  design3.chords_along_span, (loads_positive, loads_negative), spanwise_position) #testing tensile side of the design
main4(design3.moi_x_list, design3.trapezoid, design3.stringers, design3.chords_along_span, (loads_positive, loads_negative), spanwise_position, design3) #testing critical buckling case for the design
design3.graph_visualized_and_dash()

# moi_x_list, trapezoid, stringer_positions, span_positions_and_chord = main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.06808, 0.01, 0.015, 42, 4e-4)
# margin, max = main5(moi_x_list, trapezoid, span_positions_and_chord, (loads_positive, loads_negative),spanwise_position)
# main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.07702, 0.01, 0.01, 42, 4e-4)
# main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.06808, 0.01, 0.015, 42, 4e-4)
# margin, max = main5(moi_x_list, trapezoid, span_positions_and_chord, (loads_positive, loads_negative),spanwise_position)
# print(margin)